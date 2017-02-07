set -e
#Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

options="catchsegv"
for arg in "$@"
do
    options="$options '${arg}'"
done

>&2 echo "Running comand ${options}"
eval "${options}"

