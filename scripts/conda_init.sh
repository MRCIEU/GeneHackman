#!/bin/bash
set -e

__conda_setup="$('/bin/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/bin/etc/profile.d/conda.sh" ]; then
        . "/bin/etc/profile.d/conda.sh"
    else
        export PATH="/bin/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
