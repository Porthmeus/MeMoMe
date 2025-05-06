import argparse
import libsbml
import cobra
import pandas as pd
import os


def read_list_from_file(file_path):
    """Reads a list of items from a file, one item per line."""
    with open(file_path, "r") as file:
        items = [line.strip() for line in file.readlines()]
    return items


def shrink_model(
    cobra_model, desired_biomass_components, required_excretions, desired_medium, args
):
    """Simplifies the model based on given criteria."""
    # Keep a copy of the old medium and old biomass function
    old_medium = cobra_model.medium.copy()
    biomass_func = cobra_model.reactions.get_by_id(args.biomass_id)
    old_biomass = biomass_func.copy()

    # Create a reduced version of the biomass reaction
    reactions_to_add = cobra.core.DictList()
    new_react = cobra.core.reaction.Reaction("biomass_simplified")
    new_react.name = "biomass with less metabolic requirement"
    new_react.lower_bound = 0
    new_react.upper_bound = 1000
    new_react.reversibility = False
    new_react_mets = dict()
    for met in desired_biomass_components:
        new_react_met = cobra_model.metabolites.get_by_id(met)
        new_react_mets[new_react_met] = cobra_model.reactions.get_by_id(
            args.biomass_id
        ).metabolites[cobra_model.metabolites.get_by_id(met)]
    new_react.add_metabolites(new_react_mets)
    reactions_to_add.add(new_react)
    cobra_model.add_reactions(reactions_to_add)

    # Remove old biomass reaction
    cobra_model.remove_reactions(
        cobra.core.DictList([cobra_model.reactions.get_by_id(args.biomass_id)])
    )

    # Set the new biomass reaction as objective
    cobra_model.objective = "biomass_simplified"

    # Calculate optimal growth with new biomass function
    max_growth = cobra_model.slim_optimize()

    # Obtain a minimal medium able to support at least 10% of the maximal biomass growth rate
    required_medium = cobra.medium.minimal_medium(cobra_model, max_growth / 10)

    # Close all the medium's influxes
    for reac_id in cobra_model.medium.keys():
        cobra_model.reactions.get_by_id(reac_id).lower_bound = 0.0

    # Include the medium components specified by the user
    for reac_id in desired_medium.keys():
        required_medium[reac_id] = desired_medium[reac_id]

    # Set the new medium
    for reac_id in required_medium.keys():
        cobra_model.reactions.get_by_id(reac_id).lower_bound = -required_medium[reac_id]

    # Use fva to calculate the excretion and uptakes that the model can achieve while keeping 10% optimality
    fva_sol = cobra.flux_analysis.flux_variability_analysis(
        cobra_model, fraction_of_optimum=0.1
    )
    active_fva_ids = list(
        set(fva_sol[fva_sol["minimum"] < -1e-06].index).union(
            set(fva_sol[fva_sol["maximum"] > 1e-06].index)
        )
    )
    active_fva_exchanges = sorted([r_id for r_id in active_fva_ids if "EX_" in r_id])
    possible_excretions = sorted(
        [r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].maximum > 1e-06]
    )
    possible_uptakes = sorted(
        [r_id for r_id in active_fva_exchanges if fva_sol.loc[r_id].minimum < -1e-06]
    )

    # Select all the sinks and demands that can carry a flux
    boundary_ids = {r.id for r in cobra_model.boundary}
    other_reacs_like_sinks_or_demands = sorted(
        [
            r_id
            for r_id in boundary_ids.intersection(set(active_fva_ids))
            if r_id not in set(possible_excretions).union(set(possible_uptakes))
        ]
    )

    # Remove from the model all the boundary reactions that are not feasible uptakes, excretions, sinks or demands
    reacs_to_remove = list()
    for r in cobra_model.boundary:
        if r.id not in set(possible_excretions).union(set(possible_uptakes)).union(
            other_reacs_like_sinks_or_demands
        ):
            reacs_to_remove.append(r)
    cobra_model.remove_reactions(reacs_to_remove)

    # Obtain a small list of reactions that support a 10% maximal biomass production(not proven to be minimal, but with 5 iterations, and using pFBA,the total number of reactions should be small). The excretion of all the compounds required by the user will be guaranteed. All medium reactions required by the user will be kept
    for i in range(5):
        cobra_model.reactions.get_by_id("biomass_simplified").lower_bound = (
            0.1 * cobra_model.slim_optimize()
        )
        reactions_to_keep = set(required_excretions).union(set(required_medium.keys()))
        for reac in required_excretions:
            cobra_model.objective = reac
        pfba_sol = cobra.flux_analysis.pfba(cobra_model)
        pfba_sol = pfba_sol.to_frame()
        reactions_to_keep = reactions_to_keep.union(
            set(pfba_sol[abs(pfba_sol["fluxes"]) > 1e-06].index)
        )
        if pfba_sol["fluxes"][reac] <= 1e-06:
            raise Exception(
                "warning: required product " + reac + " cannot be produced by the model"
            )
        cobra_model.reactions.get_by_id("biomass_simplified").lower_bound = 0.0
        cobra_model.objective = "biomass_simplified"
        reacs_to_remove = list()
        for reac in cobra_model.reactions:
            if reac.id not in reactions_to_keep:
                reacs_to_remove.append(reac)
        cobra_model.remove_reactions(reacs_to_remove)
        # remove all unused metabolites and genes
        cobra_model, removed_metabolites = (
            cobra.manipulation.delete.prune_unused_metabolites(cobra_model)
        )
        model_reac_ids = [r.id for r in cobra_model.reactions]
        genes_to_delete = list()
        for gene in cobra_model.genes:
            useful_gene = False
        for r in gene.reactions:
            if r.id in model_reac_ids:
                useful_gene = True
        if not useful_gene:
            genes_to_delete.append(gene)
        cobra.manipulation.remove_genes(cobra_model, gene_list=genes_to_delete)

    # Set new id and name to the model
    cobra_model.id = cobra_model.id + "_shrunk"
    cobra_model.name = cobra_model.name + "after shrinking model's content"

    # This avoids problems that would occur when saving the model
    for reac in cobra_model.reactions:
        reac.lower_bound = float(reac.lower_bound)
        reac.upper_bound = float(reac.upper_bound)

    return cobra_model


def adjust_sbml_annotations(sbml_model, keep_identifiers):
    keep_identifiers_lower = [
        identifier.lower() for identifier in keep_identifiers
    ]  # Normalize to lowercase for comparison

    for species in sbml_model.getListOfSpecies():
        annotation = species.getAnnotation()

        if not annotation:
            continue

        rdf_index = -1
        for i in range(annotation.getNumChildren()):
            # print(i)
            child = annotation.getChild(i)
            # print(child)
            if child.getName() == "RDF":
                rdf_index = i
                # print(rdf_index)
                break

        if rdf_index == -1:  # RDF node not found
            continue

        rdf_node = annotation.getChild(rdf_index)

        for i in range(rdf_node.getNumChildren()):
            desc_node = rdf_node.getChild(i)
            if desc_node.getName() == "Description":
                # print("desc_node: " + str(desc_node))
                for j in range(desc_node.getNumChildren()):
                    is_node = desc_node.getChild(j)
                    # print("is_node: " + str(is_node))
                    for k in range(is_node.getNumChildren()):
                        bag_node = is_node.getChild(k)
                        # print("bag_node: " + str(is_node))
                        li_indices_to_remove = []
                        for l in range(bag_node.getNumChildren()):
                            li_node = bag_node.getChild(l)
                            li_attributes = li_node.getAttributes()
                            for ll in range(0, li_attributes.getLength()):
                                # print(f"Attribute name: {li_attributes.getName(ll)}")
                                attr_name = li_attributes.getName(ll)
                                if attr_name == "resource":
                                    resource = li_attributes.getValue(ll)
                                    # print(f"Attribute value: {resource}")
                                    if not any(
                                        identifier in resource
                                        for identifier in keep_identifiers
                                    ):
                                        # Mark index for removal if it doesn't contain any of the keep identifiers
                                        li_indices_to_remove.append(l)

                        # Remove the marked <rdf:li> elements by index, in reverse order
                        for index in reversed(li_indices_to_remove):
                            bag_node.removeChild(index)

        # Update the annotation with the modified RDF
        species.setAnnotation(annotation)


def main():
    # Setup argparse to handle input arguments

    parser = argparse.ArgumentParser(
        description="Modify a metabolic model based on specified parameters."
    )

    parser.add_argument(
        "--model_file_path",
        type=str,
        default="../tests/dat/e_coli_vmh.xml",
        help="Path to the model file.",
    )

    parser.add_argument(
        "--biomass_id",
        type=str,
        default="biomass525",
        help="ID of the biomass reaction in the model.",
    )

    parser.add_argument(
        "--desired_biomass_components_file",
        type=str,
        default="script_accessories/desired_biomass_components.txt",
        help="File containing desired biomass components.",
    )

    parser.add_argument(
        "--required_excretions_file",
        type=str,
        default="script_accessories/required_excretions.txt",
        help="File containing required excretions.",
    )

    parser.add_argument(
        "--desired_medium_file",
        type=str,
        default="script_accessories/desired_medium.csv",
        help="CSV file containing desired medium components and constraints.",
    )

    # New argument for specifying identifier types to keep
    parser.add_argument(
        "--keep_identifiers",
        nargs="+",
        required=True,
        help="List of identifier types to keep, e.g., vmh chebi bigg",
    )

    parser.add_argument(
        "--new_model_id",
        type=str,
        default="minimal_model",
        help="New ID for the modified model.",
    )

    parser.add_argument(
        "--new_model_name",
        type=str,
        default="minimal_model",
        help="New name for the modified model.",
    )

    parser.add_argument(
        "--output_file",
        type=str,
        default="../tests/dat/modified_model.xml",
        help="Path for the output modified model.",
    )

    args = parser.parse_args()

    # Load the model
    cobra_model = cobra.io.read_sbml_model(args.model_file_path)

    # Read input files
    desired_biomass_components = read_list_from_file(
        args.desired_biomass_components_file
    )
    required_excretions = read_list_from_file(args.required_excretions_file)
    desired_medium = pd.read_csv(args.desired_medium_file, header=0)
    desired_medium = pd.Series(
        [float(val) for val in desired_medium.iloc[:, 1]],
        index=desired_medium.iloc[:, 0],
    ).to_dict()

    # Call the shrink_model function
    cobra_model = shrink_model(
        cobra_model,
        desired_biomass_components,
        required_excretions,
        desired_medium,
        args,
    )

    # Create a temporary file to save the modified Cobra model as SBML
    cobra.io.write_sbml_model(
        cobra_model, args.output_file.replace(".xml", "_original_sbml.xml")
    )

    # Load the temporary SBML file with libSBML for annotation adjustment
    sbml_document = libsbml.readSBML(
        args.output_file.replace(".xml", "_original_sbml.xml")
    )
    sbml_model = sbml_document.getModel()
    print(sbml_model)

    # Adjust annotations based on the user's specified identifier types to keep
    adjust_sbml_annotations(sbml_model, args.keep_identifiers)

    # Save the final SBML model
    libsbml.writeSBMLToFile(sbml_document, args.output_file)

    os.remove(args.output_file.replace(".xml", "_original_sbml.xml"))

    print(f"Final model with adjusted annotations saved to {args.output_file}")
    pass


if __name__ == "__main__":
    # freeze_support()
    main()
