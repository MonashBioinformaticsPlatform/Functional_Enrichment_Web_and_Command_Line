library(stringr)
library(yaml)

get_script <- function(filename) {
    lines <- readLines(filename, warn = FALSE)
    lines <- str_trim(lines, "right")

    state <- "main"
    blank <- ""
    need_blank <- FALSE
    comment_code <- FALSE
    output <- c()
    show <- function(...) output[length(output) + 1] <<- paste0(...)
    for (line in lines) {
        if (state == "main" && line == "---") {
            state <- "start"
        } else if (state == "start" && str_detect(line, "^---")) {
            state <- "main"
            need_blank <- FALSE
        } else if (state == "main" && str_detect(line, "<!--omit-begin-->")) {
            state <- "omit"
        } else if (state == "omit" && str_detect(line, "<!--omit-end-->")) {
            state <- "main"
        } else if (state == "main" && str_detect(line, "^```.*include=FALSE")) {
            state <- "hidden_code"
        } else if (state == "main" && str_detect(line, "^```")) {
            state <- "code"
            comment_code <- !str_detect(line, "^``` *\\{ *[rR]")
            show("")
        } else if (state %in% c("code", "hidden_code") && str_detect(line, "^```")) {
            state <- "main"
            need_blank <- TRUE
            blank <- ""
        } else if (state == "main") {
            # Markdown links [text](url) read better in a comment as "text (url)".
            clean_line <- str_replace_all(line, "\\[([^\\]]*)\\]\\(([^)]*)\\)", "\\1 (\\2)")
            clean_line <- str_replace_all(clean_line, "`|<details>|</details>|<summary>|</summary>|\\{\\..*\\}", "")
            # Strip HTML comments and empty spacer paragraphs wherever they fall on the line.
            clean_line <- str_replace_all(clean_line, "<!--.*?-->", "")
            clean_line <- str_replace_all(clean_line, "<p>&nbsp;</p>", "")
            clean_line <- str_trim(clean_line, "right")
            if (clean_line == "" || clean_line == "\\" || clean_line == "***" || str_detect(clean_line, "^:::+$")) {
                # Blank line, or a Pandoc fenced-div marker (::: / :::: {.class}) with no
                # meaning outside the book's HTML output.
                need_blank <- TRUE
            } else if (str_detect(clean_line, "^#+")) {
                # Heading
                clean_line <- str_replace_all(clean_line, "\\{.*\\}", "")
                clean_line <- str_trim(clean_line, "right")
                depth <- nchar(str_match(clean_line, "^#*")[1, 1])
                if (depth == 1) clean_line <- paste0(clean_line, " ================")
                if (depth == 2) clean_line <- paste0(clean_line, " --------")
                if (depth <= 2) {
                    show("")
                    show("")
                } else {
                    show(blank)
                }
                show(clean_line)
                need_blank <- TRUE
                blank <- "#"
            } else {
                if (need_blank) {
                    show(blank)
                }

                indent <- nchar(str_match(clean_line, "^ *")[1, 1])
                exdent <- nchar(str_match(clean_line, "^[^A-Za-z0-9]*")[1, 1])
                wrap <- str_wrap(clean_line, indent = indent, exdent = exdent, width = 80) |> str_split_1("\n")
                for (item in wrap) show("# ", item)
                need_blank <- FALSE
                blank <- "#"
            }
        } else if (state == "code") {
            # knitr::opts_knit$set(root.dir = ...) only redirects knitr's own chunk
            # evaluation during a real render - it has no effect when this line is run
            # as plain R code, so swap it for a real setwd() that actually works here.
            line <- str_replace(line, "knitr::opts_knit\\$set\\(root\\.dir = (.*)\\)", "setwd(\\1)")
            if (comment_code) line <- paste0("# ", line)
            show(line)
        } else if (state == "hidden_code" && str_detect(line, "knitr::opts_knit\\$set\\(root\\.dir")) {
            # This line is otherwise dropped along with the rest of its include=FALSE
            # chunk, but its real-setwd() equivalent is needed so later chapters that
            # assume the project root (e.g. their own "data/R_data/..." paths)
            # resolve correctly once this chapter's chunks are done.
            show(str_replace(line, "knitr::opts_knit\\$set\\(root\\.dir = (.*)\\)", "setwd(\\1)"))
        }
    }

    output
}


# filenames <- list.files(pattern="[0-9].*\\.Rmd") |> sort() |>
#    setdiff(c("01-01-setup.Rmd", "07-Acknowledgements.Rmd", "08-session_info.Rmd"))

filenames <- read_yaml("_bookdown.yml", readLines.warn = FALSE)$rmd_files |>
    setdiff(c(
        "index.Rmd",
        "01-overview.Rmd",
        "01-02-setup.Rmd",
        "part-day1.Rmd",
        "02-recap.Rmd",
        "03-stats.Rmd",
        "04-example-dataset.Rmd",
        "05-genelists.Rmd",
        "06-1-gprofiler-web.Rmd",
        # "06-2-gprofiler2.Rmd",
        "07-1-gsea-web.Rmd",
        # "07-2-fgsea.Rmd",
        "08-1-string-web.Rmd",
        # "08-2-stringdb.Rmd",
        "09-1-reactome-web.Rmd",
        # "09-2-reactomepa.Rmd",
        "part-day2.Rmd",
        # "10-clusterprofiler.Rmd",
        # "11-novel-species-FEA.Rmd",
        "12-uncertainties.Rmd",
        "13-reporting.Rmd",
        "14-resources.Rmd"
    ))

output <- lapply(filenames, get_script) |> unlist()

writeLines(output, "docs/workshop.R")

# source("scripts/make_R_script.R")
