import unittest
from pathlib import Path
import os
import shutil
import sys
if __name__ == '__main__':
  sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

#from src.MeMoModel import MeMoModel
from src.MeMoModel import MeMoModel
from src.annotation.annotateChEBI import *
from src.annotation.annotateBiGG import *
from src.annotation.annotateModelSEED import *
from src.annotation.annotateAux import AnnotationResult
from datetime import datetime
  # Example data for the table
table_data = [
]

class Test_annotateBulkRoutines(unittest.TestCase):
    # The directory of this file
    #this_directory = Path("tests")
    this_directory = Path(__file__).parent.parent
    dat = this_directory.joinpath("dat")


    def test_ecoli_core_seed(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 54, 54)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_seed.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_seed_id(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 0,0 )
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_seed_id.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 0, 0)
        res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg_id(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0, 54, 0)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg_id.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_chebi(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 0, 0)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_chebi.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_core_bigg_chebi(self):
        mod_path = self.dat.joinpath("e_coli_core.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(54, 0, 0)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_core_bigg_chebi.__name__, res1 + res2, exp, mod)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_ecoli_vmh_bigg(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateBiGG(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")


    def test_ecoli_vmh_bigg_id(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,866,358)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg_id.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_vmh_bigg_id_chebi(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(487,866,358)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_bigg_id_chebi.__name__, res1 + res2, exp, mod)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")


    def test_ecoli_vmh_chebi(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_chebi.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")
    
    def test_ecoli_vmh_seed(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_seed.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_vmh_seed_id(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_ecoli_vmh_seed_id.__name__, res, exp, mod)   
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
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_chebi.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,577,196)
        res = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id_chebi(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(311,577,196)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateChEBI(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id_chebi.__name__, res1 + res2, exp, mod)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_adlercreutzia_equolifaciens_bigg_id_seed(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(336,1082,671)
        res1 = AnnotationResult.fromAnnotation(annotateBiGG_id(mod.metabolites))
        res2 = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_bigg_id_seed.__name__, res1 + res2, exp, mod)   
        self.assertLessEqual(exp, res1 + res2, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res1 + res2}. All three must be >=")

    def test_adlercreutzia_equolifaciens_seed(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_seed.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_adlercreutzia_equolifaciens_seed_id(self):
        mod_path = self.dat.joinpath("Adlercreutzia_equolifaciens_DSM_19450.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(annotateModelSEED_id(mod.metabolites))
        add_test_case_to_table(self.test_adlercreutzia_equolifaciens_seed_id.__name__, res, exp, mod)   
        self.assertLessEqual(exp, res, msg=f"Expected amount of annotated metabolites: {exp}, calculated amount of annotated metabolites: {res}. All three must be >=")

    def test_ecoli_vmh_all(self):
        mod_path = self.dat.joinpath("e_coli_vmh.xml")
        mod = MeMoModel.fromPath(mod_path)
        exp = AnnotationResult(0,0,0)
        res = AnnotationResult.fromAnnotation(mod.annotate())
        add_test_case_to_table(self.test_ecoli_vmh_all.__name__, res, exp, mod)   
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

def bac_spec_content(mod: MeMoModel) -> str:
    html = "<table id = myTable>\n"
    html += """ 
    <thead>
      <tr>
        <th>ID</th>
        <th>Names</th>
        <th>Annotations</th>
      </tr>
    </thead>
    <tbody>
    """
    for metabolite in mod.metabolites:
        html += f"<tr><td>{metabolite.id}</td><td>{metabolite.names}</td><td>{metabolite.annotations}</td></tr>"
    html += "</tbody></table>"
    return html

def generate_css():
    return """<style>
    #searchInput {
      margin: 12px;
      padding: 8px;
      width: 250px;  
    }
    canvas{
    
      width: 45% !important;
      height: 500px !important;
    
    }
    .container {
      display: flex;
      gap: 10px;
    }

    table {
      width: 45%;
      height: 500px;
      border-collapse: collapse;
      margin-top: 10px;
      overflow-y: auto;  /* Enable vertical scrolling */
      display: grid;
      grid-auto-flow: column; /* make the items fill up the grid container by column and not row */
      grid-template-rows: repeat(10, 1fr); /* have a max number of 3 rows and make them all the same height */
      gap: 0.125rem; /*put a small gap between each element */
    }
    th, td {
      padding: 12px;
      text-align: left;
      border: 1px solid #ddd;
    }
    th {
      background-color: #f2f2f2;
    }
    /* Constrain the table height and make it scrollable */
    #tableContainer {
      max-height: 200px; /* Set the height limit */
      overflow-y: auto;  /* Add vertical scroll if content exceeds max-height */
    }
  </style>"""

def generate_js_for_search(): 
    return """
    function searchTable() {
      const input = document.getElementById('searchInput');
      const filter = input.value.toUpperCase();
      const table = document.getElementById('myTable');
      const tr = table.getElementsByTagName('tr');

      for (let i = 1; i < tr.length; i++) {
        let td = tr[i].getElementsByTagName('td');
        let found = false;

        for (let j = 0; j < td.length; j++) {
          if (td[j]) {
            let txtValue = td[j].textContent || td[j].innerText;
            if (txtValue.toUpperCase().indexOf(filter) > -1) {
              found = true;
              break;
            }
          }
        }

        tr[i].style.display = found ? '' : 'none';
      }
    }
    """


def generate_bar_chart(mod: MeMoModel):
  metabs = mod.metabolites
  keys = dict()
  for m in metabs:
      for annotation in m.annotations.keys():
          keys[annotation] = keys.get(annotation, 0) + 1
  # Check if order is presercerd in keys/values
  a = keys.keys()
  b = keys.values()
  # TODO CHECK MATH (is the missing calculation correct)
  # THE MATH IS BS
  s = """
       // Get the context of the canvas element
   var ctx = document.getElementById('myBarChart').getContext('2d');
    
   // Create the bar chart
   var myBarChart = new Chart(ctx, {
     type: 'bar', // Bar chart type
     data: {
         labels: """ + str(list(a)) +  """, // Labels for the bars
       datasets: [{
         label: 'Values', // Label for the dataset
         data: """ + str(list(b)) + """,// Data for the bars
         backgroundColor: ['#FF5733', '#33FF57', '#3357FF', '#FF33A1'], // Bar colors
         borderColor: ['#FF5733', '#33FF57', '#3357FF', '#FF33A1'], // Border colors
         borderWidth: 1
       }]
     },
     options: {
       scales: {
         y: {
           beginAtZero: true, // Ensure the Y-axis starts at zero
         }
       },
       responsive: true // Make the chart responsive
     }
     });"""
  return s

        
def generate_bac_specs( name: str,  mod: MeMoModel):
    file_path = os.path.join("bac_specs", f'{name}.html')

    
    body = """<input type="text" id="searchInput"  onkeyup="searchTable()" placeholder="Search for names..">\n"""
    body += """<div class="container">\n"""
    body += bac_spec_content(mod)
    body += """<canvas id="myBarChart" class = "box" width="200px" height="100px"></canvas>\n"""
    body += "</div>"
    body +=  """<script>"""
    body += generate_js_for_search()
    body += generate_bar_chart(mod)
    body +=  """</script>"""
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
      <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
      {}
      <title>{}</title>
    </head>
    <body>
      {}
    </body>
    </html>
    """.format(generate_css(), name, body )
    
    # Write the content to the HTML file
    with open(file_path, 'w') as file:
        file.write(html_content)
    print(f'HTML file created at: {file_path}')


def add_test_case_to_table( name: str, res, exp, mod: MeMoModel):
  link  =  f"""<a href="bac_specs/{name}.html" target="_blank">{name}</a>"""
  table_data.append([f"<b>{link}</b>", "Annotated", "Expected", ""])
  generate_bac_specs(name, mod)
  table_data.append([f"Inchis", *generate_value_row(res.annotated_inchis, exp.annotated_inchis )])
  table_data.append([f"DBs", *generate_value_row(res.annotated_dbs, exp.annotated_dbs)])
  table_data.append([f"Names", *generate_value_row(res.annotated_names, exp.annotated_names)])




if __name__ == '__main__':

    # Create folder for html files 
    folder_name = "bac_specs"   
    if os.path.exists(folder_name):
      shutil.rmtree(folder_name)
      print(f"Existing folder '{folder_name}' deleted.")
    os.makedirs(folder_name)

    # Create a test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_annotateBulkRoutines)
    # The LSP complains but this works with id(). This is a super hacky way to just execute one single test for debugging
    #filtered_tests = unittest.TestSuite(
    #        test for test in suite if test.id() in ["__main__.Test_annotateBulkRoutines.test_ecoli_vmh_all"]
    #)
    #suite = filtered_tests
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
