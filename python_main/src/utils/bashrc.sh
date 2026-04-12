# ----------------------------- MATLAB Settings -----------------------------

## add to path
export PATH="/Applications/MATLAB_R2025b.app/bin:$PATH"

# MatLab cli repl
matlabrepl(){
    matlab -nodisplay 
}

## Matlab source function
# srcmatlab(){
#     matlab -nodesktop -nosplash -noFigureWindows -batch  "${@%.m}"  
# }
srcmatlab(){
    local abs                                                                                                    
    abs=$(realpath "$1")
    local dir
    dir=$(dirname "$abs")
    local base                                                                                                   
    base=$(basename "${abs%.m}")
    matlab -nodesktop -nosplash -noFigureWindows -batch "addpath('$dir'); $base"                                 
}               
