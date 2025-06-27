def create_deletion_script():
    with open('delete_generated_files.sh', 'w') as script:
        script.write("#!/bin/bash\n")
        script.write("rm *_unmodified_peptides_usi.txt\n")
        script.write("rm *_modified_peptides_usi.txt\n")


create_deletion_script()

print("files deleted successfully")