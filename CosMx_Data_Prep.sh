# Note that this code works for a single slide. If exporting
# multiple slides, some modification might be necessary.

# cd to parent directory

# create a directory adjacent to RawData and flatFiles
mkdir -p sample_dir_formatted && cd $_

# Add flat files
for file in $(ls ../flatFiles/*/*csv*); do cp $file ./; done

# Add folders
cp -r ../RawFiles/*/*/CellStatsDir/CellComposite ./
cp -r ../RawFiles/*/*/CellStatsDir/CellOverlay ./

mkdir -p CellLabels
for file in $(ls ../RawFiles/*/*/CellStatsDir/FOV*/CellLabels*); do cp $file ./CellLabels/ ; done
rm ./CellLabels/._Cell*
  
mkdir -p CompartmentLabels
for file in $(ls ../RawFiles/*/*/CellStatsDir/FOV*/CompartmentLabels*); do cp $file ./CompartmentLabels/ ; done