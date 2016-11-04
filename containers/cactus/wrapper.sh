set -e
#Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

options=""
for arg in "${@:2}"
do
    options="$options '${arg}'"
done

>&2 echo "Running comand " /home/cactus/bin/$1 "${options}"
eval /home/cactus/bin/$1 "${options}"

