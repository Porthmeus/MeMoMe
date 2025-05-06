# Porthmeus
# 30.06.24

import unittest
import cobra as cb
from src.removeDuplicateMetabolites import *


class Test_removeDuplicateMetabolites(unittest.TestCase):

    def createToyModel(self):

        toy = cb.Model("toy")
        reactionAB = cb.Reaction("R_AB")
        reactionAB.name = "A+B= C"
        reactionAB.subsystem = "test"
        reactionAB.lower_bound = -30.0
        reactionAB.upper_bound = 40.0

        reactionAD = cb.Reaction("R_AD")
        reactionAD.name = "A+D=C"
        reactionAD.subsystem = "test"
        reactionAD.lower_bound = 0.0  # This is the default
        reactionAD.upper_bound = 1000.0  # This is the default

        reactionDE = cb.Reaction("R_DE")
        reactionDE.name = "D+E = F"
        reactionDE.subsystem = "test"
        reactionDE.lower_bound = 0.0  # This is the default
        reactionDE.upper_bound = 1000.0  # This is the default

        A = cb.Metabolite("A", formula="A", name="A", compartment="c")
        B = cb.Metabolite("B", formula="B", name="B", compartment="c")
        C = cb.Metabolite("C", formula="C", name="C", compartment="c")
        D = cb.Metabolite("D", formula="D", name="D", compartment="c")
        E = cb.Metabolite("E", formula="E", name="E", compartment="c")
        F = cb.Metabolite("F", formula="F", name="F", compartment="c")

        reactionAB.add_metabolites({A: -1.0, B: -1.0, C: 1.0})
        reactionAD.add_metabolites({A: -1.0, D: -1.0, C: 1.0})
        reactionDE.add_metabolites({E: -1.0, D: -1.0, F: 1.0})

        toy.add_reactions([reactionAB, reactionAD, reactionDE])
        reactionAB.gene_reaction_rule = "g1 and g2"
        reactionAD.gene_reaction_rule = "g1 and g3"
        reactionDE.gene_reaction_rule = "g4 and g5"
        # add gene names
        for gene in toy.genes:
            gene.name = gene.id.replace("g", "gene")
        return toy

    def test_findCommonReaction(self):
        toy = self.createToyModel()
        met1 = toy.metabolites.get_by_id("B")
        met2 = toy.metabolites.get_by_id("D")
        # find common reactions - ignore if the fact that one might be reversible while the other is not
        rxns = findCommonReactions(met1, met2)
        self.assertTrue(len(rxns) == 1)
        self.assertTrue(rxns[0][0].id == "R_AB")
        self.assertTrue(rxns[0][1].id == "R_AD")

        # find only reactions if they are really the same, including the reversible information
        rxns = findCommonReactions(met1, met2, reversible_is_same=False)
        self.assertTrue(len(rxns) == 0)

    def test_mergeDuplicateReactions(self):
        """Test merging of duplicate reactions in a model"""
        toy = self.createToyModel()
        rxn1 = toy.reactions.get_by_id("R_AB")
        rxn2 = toy.reactions.get_by_id("R_AD")
        no_rxns_pre = len(toy.reactions)
        mergeReactions(rxn1, rxn2)
        self.assertTrue(no_rxns_pre != len(toy.reactions))
        self.assertTrue(len(toy.reactions) == 2)

    def test_exchangeMetabolite(self):
        toy = self.createToyModel()
        met1 = toy.metabolites.get_by_id("B")
        met2 = toy.metabolites.get_by_id("D")
        rxns_forms = [x.reaction for x in toy.reactions]
        no_mets = len(toy.metabolites)
        res = exchangeMetabolite(met1, met2, prune=False)
        # print(res)
        self.assertTrue(
            all(
                [
                    x == y
                    for x, y in zip(rxns_forms, [z.reaction for z in toy.reactions])
                ]
            )
            == False
        )
        # check if the metabolite is still associated to the model
        self.assertTrue(len(toy.metabolites) == no_mets)
        # check if pruning workd
        exchangeMetabolite(met1, met2, prune=True)
        self.assertTrue(len(toy.metabolites) == no_mets - 1)

    def test_mergeAndRemove(self):
        toy = self.createToyModel()
        met1 = toy.metabolites.get_by_id("B")
        met2 = toy.metabolites.get_by_id("D")
        res = mergeAndRemove(met1, met2)
        # print(res)
        self.assertTrue(len(toy.reactions) == 2)
        self.assertTrue(len(toy.metabolites) == 5)
