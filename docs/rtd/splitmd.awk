BEGIN {
    found_header = 0
    num_lines = 0
    first_header = 1
    num_headers = 0
}

/^=+$/ {
    found_header = 1
    num_headers += 1
    if (first_header == 1) {
        first_header = 0
    } else if (first_header == 0 && found_header) {
        print_section(section_lines, num_lines - 1, headers, num_headers)
        section_lines[1] = section_lines[num_lines]
        num_lines = 1
    }
}

{
    num_lines += 1
    section_lines[num_lines] = $0
}

function print_section(lines, n_lines, h, n_headers) {
    filename = tolower(lines[1])
    gsub(/ /, "_", filename)
    gsub(/'/, "", filename)
    h[n_headers] = filename
    filename = filename ".rst"
    pipe = "pandoc -f markdown -t rst -o " filename
    for (i = 1; i < n_lines; i++) {
        print lines[i] | pipe
    }
    close(pipe)
}

END {
    print_section(section_lines, num_lines + 1, headers, num_headers)

    # generate index.rst
    print ".. toctree::" > "index.rst"
    print "   :maxdepth 2" > "index.rst"

    for (i = 1; i <= num_headers; i++) {
        print "   " headers[i] > "index.rst"
    }

    close("index.rst")
}
