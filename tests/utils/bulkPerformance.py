import unittest
from pathlib import Path
import cobra as cb
import pandas as pd

#from src.MeMoModel import MeMoModel
from src.MeMoMetabolite import MeMoMetabolite
from src.MeMoModel import MeMoModel
from src.annotateChEBI import *
from src.annotateBiGG import *
from src.annotateModelSEED import *
from src.annotateAux import AnnotationResult
from datetime import datetime

  # Example data for the table
table_data = [
]

class Test_annotateBulkRoutines(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent
    dat = this_directory.joinpath("dat")


    def test_ecoli_core_seed(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 54, 54)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_seed.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_seed_id(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 0,0 )
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_seed_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 0, 0)
        res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg_id(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 54, 0)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_chebi(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 0, 0)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_chebi.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg_chebi(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 0, 0)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg_chebi.__name__, res1 + res2, exp)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_ecoli_vmh_bigg(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(999,999,999)
        res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")


    def test_ecoli_vmh_bigg_id(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,866,358)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_vmh_bigg_id_chebi(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,866,358)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg_id_chebi.__name__, res1 + res2, exp)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")


    def test_ecoli_vmh_chebi(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(999,999,999)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_chebi.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")
    
    def test_ecoli_vmh_seed(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_seed.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_vmh_seed_id(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_seed_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    # ############################ AT THE MOMENT NOTHING GETS ANNOTATED FOR RECON ANYWAY SO WE IGNORE IT FOR NOW ####################################33
    #def test_recon3D_301_bigg(self):
    #    # load the e.coli core model and bulk annotate the metabolites. Check if any annoation tkes place (Chebi should cover all metabolites)
    #    mod_path = self.dat.joinpath("Recon3DModel_301.xml")
    #    mod = MeMoModel.fromPath(mod_path)
    #    exp = AnnotationResult(0, 54, 0)
    #    res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
    #    add_test_case_to_table(self.test_recon3D_301_bigg.__name__, res, exp)   
    #    self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    #def test_recon3D_301_bigg_id(self):
    #    # load the e.coli core model and bulk annotate the metabolites. Check if any annoation tkes place (Chebi should cover all metabolites)
    #    mod_path = self.dat.joinpath("Recon3DModel_301.xml")
    #    mod = MeMoModel.fromPath(mod_path)
    #    exp = AnnotationResult(999,999,999)
    #    res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
    #    add_test_case_to_table(self.test_recon3D_301_bigg_id.__name__, res, exp)   
    #    self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    #def test_recon3D_301_chebi(self):
    #    # load the e.coli core model and bulk annotate the metabolites. Check if any annoation tkes place (Chebi should cover all metabolites)
    #    mod_path = self.dat.joinpath("Recon3DModel_301.xml")
    #    mod = MeMoModel.fromPath(mod_path)
    #    exp = AnnotationResult(999,999,999)
    #    res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
    #    add_test_case_to_table(self.test_recon3D_301_chebi.__name__, res, exp)   
    #    self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_chebi(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_chebi.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(999,999,999)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,577,196)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id_chebi(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,577,196)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id_chebi.__name__, res1 + res2, exp)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id_seed(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,577,196)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id_seed.__name__, res1 + res2, exp)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_adlercreutzia_equolifaciens_seed(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_seed.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_seed_id(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_seed_id.__name__, res, exp)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

def generate_html_table(data):
    """
    Generate HTML table from a list of lists (2D array) representing the table data.
    
    Parameters:
    data (list of lists): 2D array representing the table data.

    Returns:
    str: HTML code for the table.
    """
    red_string = "style=\"color: white; background-color: red;\""
    green_string = "style=\"color: white; background-color: green;\""

    html = "<table>\n"
    for row in data:
        print(row)
        html += "  <tr>\n"
        for cell in row[0:3]:
            if row[3] == "red":
              html += f"    <td {red_string}>{cell}</td>\n"
            elif row[3] == "green" :
              html += f"    <td {green_string}>{cell}</td>\n"
            else:
              html += f"    <td>{cell}</td>\n"

        html += "  </tr>\n"
    html += "</table>"
    return html


def save_html_table(html_table, filename_prefix='performace_table'):
    """
    Save HTML table to an HTML file with a filename containing the current date and time.
    
    Parameters:
    html_table (str): HTML code for the table.
    filename_prefix (str): Prefix for the filename (default is 'table').
    """
    # Get current date and time
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    
    # Construct filename
    filename = f"{filename_prefix}.html"#_{timestamp}.html"
    
    # Save HTML table to file
    with open(filename, 'w') as file:
        file.write(html_table)
# Print the generated HTML table

def generate_value_row(res_value: int, exp_value: int) -> tuple[str, str, str]:
  """
  Generate a row for displaying result and expected values with color indication.

  Parameters:
  - res_value (int): The actual result value.
  - exp_value (int): The expected value.

  Returns:
  tuple[str, str, str]: A tuple containing three strings:
      - The result value.
      - The expected value.
      - A string indicating the color for styling:
          - "red" if the result value is less than the expected value.
          - "green" if the result value is greater than or equal to the expected value.
  """
  if res_value < exp_value:
   return (str(res_value), str(exp_value), "red")
  else:
   return (str(res_value), str(exp_value), "green")

def add_test_case_to_table( name: str, res, exp):
  table_data.append([f"<b>{name}</b>", "Annotated", "Expected", ""])
  table_data.append([f"Inchis", *generate_value_row(res.annotated_inchis, exp.annotated_inchis )])
  table_data.append([f"DBs", *generate_value_row(res.annotated_dbs, exp.annotated_dbs)])
  table_data.append([f"Names", *generate_value_row(res.annotated_names, exp.annotated_names)])




if __name__ == '__main__':
    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_annotateBulkRoutines)
    # Create a test runner
    runner = unittest.TextTestRunner()
    
    # Run the tests and get the result
    result = runner.run(suite)


  
    # Generate HTML table
    html_table = generate_html_table(table_data)
    save_html_table(html_table)
    
    # Access the test results
    print("Number of tests run:", result.testsRun)
    print("Number of failures:", len(result.failures))
    print("Number of errors:", len(result.errors))
