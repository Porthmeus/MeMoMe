import cobra
from src.MeMoModel import *

def mergeModels(model1, model2, matches):
    matches = matches.loc[matches["Name_score"]==1.0]
    print(matches)
    print()
    print(model1.cobra_model.compartments)
    print()
    print(model2.cobra_model.compartments)
    print()
    # create an empty model
    merged_model = MeMoModel(cobra_model=cobra.Model(id_or_model="merged_model", name="merged_model"))


    pass


if __name__ == '__main__':

    tinyModel1 = "../tests/dat/my_model_with_annotations.xml"
    tinyModel2 = "../tests/dat/contiguous_pathway_model.xml"

    model1 = MeMoModel.fromPath(Path(tinyModel1))
    model2 = MeMoModel.fromPath(Path(tinyModel2))
    matches = model1.match(model2)
    mergeModels(model1, model2, matches)