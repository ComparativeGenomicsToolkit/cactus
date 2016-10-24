set -e
#Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

>&2 echo "Running comand " /home/cactus/bin/${1} ${@:2}
>&2 ls -la
exec /home/cactus/bin/${1} ${@:2}

