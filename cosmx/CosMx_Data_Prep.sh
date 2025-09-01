# Note that this code works for a single slide.
# The name of the slide needs to be the first input
# To run the script, cd to the directory containing RawFiles and flatFiles directories 

# create a directory adjacent to RawData and flatFiles
mkdir -p sample_dir_formatted/$1 && cd $_

# Add flat files
for file in $(ls ../../flatFiles/$1/*csv*); do cp $file ./; done

# Copy folders
cp -r ../../RawFiles/$1/CellComposite ./
cp -r ../../RawFiles/$1/CellOverlay ./

# Extract CellLabels from each FOV folder
mkdir -p CellLabels
for file in $(ls ../../RawFiles/$1/FOV*/CellLabels*); do cp $file ./CellLabels/ ; done

# Extract CompartmentLabels from each FOV folder 
mkdir -p CompartmentLabels
for file in $(ls ../../RawFiles/$1/FOV*/CompartmentLabels*); do cp $file ./CompartmentLabels/ ; done