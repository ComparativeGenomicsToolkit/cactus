set -e
# We forward any termination signals we receive to the underlying
# process tree. -1 is always our process group ID when running in a
# container.
sigint() {
    kill -SIGINT -- -1
}

trap sigint SIGINT

sigterm() {
    kill -SIGTERM -- -1
}

trap sigterm SIGTERM

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
eval "${options}" <&0 &
wait $!
exit $?
