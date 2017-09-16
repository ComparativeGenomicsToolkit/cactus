set -e
# Set monitor mode: give child processes a new PGID
set -m

# Find pgid of the child process
# Credit: https://stackoverflow.com/a/36820679
pgid_from_pid() {
    ps -o pgid= "$pid" 2>/dev/null | egrep -o "[0-9]+"
}

# Forward signals to the child process tree
sigint() {
    kill -SIGINT -- -$(pgid_from_pid)
    wait $pid
}

trap sigint SIGINT

sigterm() {
    kill -SIGTERM -- -$(pgid_from_pid)
    wait $pid
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
pid=$!
wait $pid
exit $?
