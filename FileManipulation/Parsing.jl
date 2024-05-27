# Looping through a file, finding expression and replacing them

function parse_file()

    # Set data file path
    file_path_old = "./FileManipulation/Literature.bib"
    file_path_new = "./FileManipulation/Literature_new.bib"
    
    # What to search for
    target_str    = "title"     # Line that contains this string should be replaced

    # Call function
    find_replace_occurence(file_path_old, file_path_new, target_str)

    return nothing
end

function find_occurence(file, target_str)
    idx = Tuple{Int,Int}[]
    for (i, line) in enumerate(eachline(file))
        if (m = match(Regex("$(target_str)"), line); m !== nothing)
            println(new_line)
            push!(idx, (i, m.offset))
        end
    end

    return idx
end

function find_replace_occurence(file_old, file_new, target_str)
    open(file_new, "w") do fnew
        for (i, line_old) in enumerate(eachline(file_old))
            if (m = match(Regex("$(target_str)"), line_old); m !== nothing)
                title_str_old = match(r"\{(.*?)\}", line_old)                                       # Replace everything inside curly braces of current line ...
                line_new = replace(line_old, "$(title_str_old[1])" =>"{$(title_str_old[1])}")       # ... with the same content enclosed in a second pair of curly braces
                println(fnew, line_new)
            else
                println(fnew, line_old)
            end
        end
    end

    # Return
    return nothing
end

# Run
parse_file()
