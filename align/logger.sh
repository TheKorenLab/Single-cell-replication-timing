log() {
    printf '%5s' ;
    echo "$(date)" ;
    echo "$@" ;
    printf '%50s\n' | tr ' ' - ;
} >&2
