#!/bin/bash
#
# Compile the BEst standalone application with MATLAB
# 
# Last revision: July, 2019
# Maintainer: Obai Bin Ka'b Ali @aliobaibk
# License: In the app folder or check GNU GPL-3.0.



############## Function to check if help was requested
function help_requested() {
    [[ "$@" =~ (^|[[:space:]])"-h"($|[[:space:]]) ]] || \
    [[ "$@" =~ (^|[[:space:]])"--help"($|[[:space:]]) ]]
}



############## Function to check a download status
function check_download_status() {
    exit_code="$1"
    file_desc="$2"
    if [[ "$exit_code" != 0 ]]; then
        echo -e "\n\n\n     ***** Failed to download [""$file_desc""]"
        echo -e "\n     - Please check (or report a bug at): https://github.com/multifunkim/best-cbrain"
        echo -e "\n     - Closing the program\n"
        exit 1
    fi
}



# Inputs
this_loc="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if [ ! -f "$this_loc"/LICENSE ]; then
    echo -e "\n     - The license file is missing, downloading now..."

    wget -q -P "$this_loc" "https://raw.githubusercontent.com/multifunkim/best-cbrain/master/LICENSE"
    check_download_status $? "license file"
fi

if [ ! -f "$this_loc"/best_cargparse.py ]; then
    echo -e "\n     - The file 'best_cargparse.py' is missing, downloading now..."

    wget -q -P "$this_loc" "https://raw.githubusercontent.com/multifunkim/best-cbrain/master/for_build/best_cargparse.py"
    check_download_status $? "best_cargparse.py"
    chmod u+x "$this_loc"/best_cargparse.py
fi

oargs="$("$this_loc"/best_cargparse.py "$@")"
if [[ $? != 0 ]]; then
    echo -e "\n\n\n     ***** Something went wrong"
    echo -e "\n     Closing the program\n"
    exit 1
  
elif help_requested "$@"; then
    echo "$oargs"
    exit 0
  
fi

output_dir="$(echo "$oargs" | cut -d '"' -f 2)"
output_name="$(echo "$oargs" | cut -d '"' -f 4)"
matlab_cmd="$(echo "$oargs" | cut -d '"' -f 6)"



# Temporary directory
tmp_dir="$(mkdir -p "$output_dir" && mktemp -d --tmpdir="$output_dir" tmp-XXXXX)"



# URLs
url_best_lib="https://api.github.com/repos/multifunkim/best-brainstorm/releases/latest"
url_bst="https://raw.githubusercontent.com/brainstorm-tools/brainstorm3/master/external/brainentropy"
url_main="https://raw.githubusercontent.com/multifunkim/best-cbrain/master/for_build/app_extra/BEst.m"



############## Cleaning function
clean () {
    echo -e "\n\n\n     ***** Doing some cleaning, PLEASE WAIT"

    rm -rf "$tmp_dir" >/dev/null 2>&1

    echo -e "\n     - All cleaning done, BYE"
}
trap clean EXIT



# BEst library directory and main function
best_lib="$tmp_dir"/lib
best_main="$tmp_dir"/BEst.m



# Download best-brainstorm and unzip
echo -e "\n     - Downloading the BEst library..."

zip_url="$(curl -s "$url_best_lib" | grep zipball_url | cut -d '"' -f 4)"
curl -s -L -o "$tmp_dir"/best-lib.zip "$zip_url"
check_download_status $? "BEst library"
unzip -q "$tmp_dir"/best-lib.zip -d "$best_lib"
mv "$best_lib"/multifunkim-best-brainstorm*/* "$best_lib"
rm -rf "$tmp_dir"/best-lib.zip "$best_lib"/multifunkim-best-brainstorm*



# Download Brainstorm program files
echo -e "\n     - Downloading the Brainstorm program files..."

wget -q -P "$best_lib"/bst-files "$url_bst"/be_install.m
check_download_status $? "be_install.m"
wget -q -P "$best_lib"/bst-files "$url_bst"/be_main.m
check_download_status $? "be_main.m"
wget -q -P "$best_lib"/bst-files "$url_bst"/panel_brainentropy.m
check_download_status $? "panel_brainentropy.m"



# Download BEst main file
echo -e "\n     - Downloading BEst main file..."

wget -q -O "$best_main" "$url_main"
check_download_status $? "BEst main file"



# Compile
echo -e "\n     - Compiling..."

app_zip="$output_name".zip
app_dir="$tmp_dir"/"${app_zip%'.zip'}"
mkdir "$app_dir"

cp "$this_loc"/LICENSE "$app_dir"

mcc_cmd='mcc -m '"$best_main"' -a '"$best_lib"' -d '"$app_dir"' -o best_samapp -R -nodisplay -v'
"$matlab_cmd" -nodisplay -nosplash -r \
'try, '"$mcc_cmd"', catch e, disp(getReport(e)), exit(1), end, exit(0);'

if [[ $? != 0 ]]; then
    echo -e "\n     - ERROR, compilation failed"
    exit 1
elif [ -z "$(ls -A "$app_dir")" ]; then
    echo -e "\n     - ERROR, empty directory:\n""$app_dir""\n"
    exit 1
else
    pushd "$app_dir" >/dev/null 2>&1
    zip -qr "$output_dir"/"$app_zip" .
    popd >/dev/null 2>&1
    echo -e "\n     - All done, check:\n""$output_dir"/"$app_zip"
fi
