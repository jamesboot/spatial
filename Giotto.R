# Script for trialling Giotto package for spatial transcriptomics
setwd('/Volumes/My Shared Files/Home/Documents/projects/GC-IM-10692/')

# Install  dependencies
#install.packages('terra')
#install.packages('pak')

# Install Giotto
#pak::pkg_install('drieslab/Giotto')

# Load Giotto
library(Giotto)


reticulate::conda_install(
  'r-reticulate',
  packages = c(
    "pandas==1.5.1",
    "networkx==2.8.8",
    "python-igraph==0.10.2",
    "leidenalg==0.9.0",
    "python-louvain==0.16",
    "python.app==1.4",
    "scikit-learn"
  )
)

reticulate::use_condaenv('r-reticulate')

reticulate::conda_list()

dir.create('Giotto_outs')
results_dir = 'Giotto_outs'  
my_python_path = '/opt/anaconda3/envs/r-reticulate/bin/python'

instrs = createGiottoInstructions(save_dir = results_dir,
                                  python_path = my_python_path)


# Create python env for Giotto
installGiottoEnvironment(force_miniconda = T, force_environment = T)

.rs.restartR()


# Creating Giotto Instructions without specifying a Python path will make 
# reticulate activate the default Giotto environment. 
default_instrs <- createGiottoInstructions()

# Extract python path information
default_python_path <- default_instrs$python_path

# Direct reticulate to use Python within the Giotto Environment
reticulate::use_python(default_python_path)

# Create python env for Giotto
installGiottoEnvironment()

reticulate::
