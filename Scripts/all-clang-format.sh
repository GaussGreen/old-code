#!/bin/bash

cd ..

MAIN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Variable that will hold the name of the clang-format command
FMT=""

# Some distros just call it clang-format. Others (e.g. Ubuntu) are insistent
# that the version number be part of the command. We prefer clang-format if
# that's present, otherwise we work backwards from highest version to lowest
# version.
for clangfmt in clang-format{,-{4,3}.{11,10,9,8,7,6,5,4,3,2,1,0}}; do
    if which "$clangfmt" &>/dev/null; then
        FMT="$clangfmt"
        break
    fi
done

# Check if we found a working clang-format
if [ -z "$FMT" ]; then
    echo "failed to find clang-format"
    exit 1
fi

# Run clang-format -i on all of the things
files_added=`git diff --cached --name-only | egrep '\.c$|\.cxx$|\.cpp$|\.cu$|\.h$|\.hxx$|\.h.in$|\.hpp$|\.java$' `
files_existent=`git diff --name-only | egrep '\.c$|\.cxx$|\.cpp$|\.cu$|\.h$|\.hxx$|\.h.in$|\.hpp$|\.java$' `
"${FMT}" -verbose -i $files_existent