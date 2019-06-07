{
    printf($1)
    for (i=2; i<=NF; i++) {
        if (i <= 6) {
            printf("\t%s", $i)
        } else {
            split($i, arr, ",")
            printf("\t%s", arr[field])
        }
    }
    printf("\n")
}
