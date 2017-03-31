set -e
options="catchsegv"
for arg in "$@"
do
    options="$options '${arg}'"
done

>&2 echo "Running comand ${options}"
eval "${options}"

