import json
import re

#example JSP file:
'''
<script language="javascript" type="text/javascript">
        var table = null;
        var dataset_files = {"compliance_threshold":"Jul. 2, 2024, 1:59 PM","total_rows":1160,"inactive_threshold":"Nov. 2, 2023, 1:59 PM","source":"ProteoSAFe","row_data":[{"last_used_millis":1730583002856,"size_MB":0.10456085205078125,"descendant_files":0,"size":109640,"attachment":"","file_descriptor":"f.MSV000096122\/ccms_parameters\/params.xml","last_used":"Nov. 2, 2024, 1:30 PM","name":"params.xml","collection":"ccms_parameters","relative_path":"ccms_parameters","type":"File"},{"last_used_millis":1729574483883,"size_MB":819.7329912185669,"descendant_files":0,"size":859552341,"attachment":"","file_descriptor":"f.MSV000096122\/ccms_peak\/RAW\/ScltlMsclSet10BRPhsFr1.mzML","last_used":"Oct. 21, 2024, 9:21 PM","name":"ScltlMsclSet10BRPhsFr1.mzML","collection":"ccms_peak","relative_path":"ccms_peak\/RAW","type":"File"},{"last_used_millis":1729574502014,"size_MB":767.3281908035278,"descendant_files":0,"size":804601925,"attachment":"",

'''
# Read the JSP file
with open('dataset_files.jsp', 'r') as file:
    jsp_content = file.read()

# Attempt to match the dataset_files variable with a more flexible regex pattern
json_match = re.search(r'dataset_files\s*=\s*({.*?})\s*;\s*', jsp_content, re.DOTALL)

# Check if the match was successful
if json_match:
    # Extract and clean up JSON-like content
    json_content = json_match.group(1)
    json_content = re.sub(r',\s*}', '}', json_content)  # Remove trailing commas before closing braces
    json_content = re.sub(r',\s*\]', ']', json_content)  # Remove trailing commas before closing brackets
    
    try:
        # Parse the JSON content
        data = json.loads(json_content)
        
        # Extract all .mzML file names
        mzml_files = [item['name'] for item in data['row_data'] if item['name'].endswith('.mzML')]
        
        # Write the filenames to a text file
        with open('mzml_files.txt', 'w') as output_file:
            for filename in mzml_files:
                output_file.write(filename + '\n')
        
        print("mzML filenames successfully written to mzml_files.txt")
    except json.JSONDecodeError as e:
        print("JSON decoding error:", e)
else:
    # Print the JSP content to help with debugging
    print("No JSON data found in JSP content. Check the structure of the file.")
