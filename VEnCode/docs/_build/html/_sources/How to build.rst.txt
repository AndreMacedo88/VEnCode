in terminal::

$ sphinx-build -b html . .\build

explanation: sphinx-build -b method source_path_of_config output_folder_name

or better yet::

$ make html

explanation: if there is already a makefile in the folder, it's quicker. You can also do::

$ make clean

to clear the build folder first.
