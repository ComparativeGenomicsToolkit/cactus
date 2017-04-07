set -e
options="catchsegv"
for arg in "$@"
do
    if [ "$arg" == "cactus-redirect" ]; then
        options="$options >"
    else
        options="$options '${arg}'"
    fi
done

>&2 echo "Running command ${options}"
eval "${options}"
