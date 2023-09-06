# Identify all source code files in (sub)folders
files <- list.files(pattern='[.](R|rmd)$', all.files=T, recursive=T,
                    full.names = T, ignore.case=T)

# Extract the source code
code = unlist(sapply(files, scan, what = 'character', quiet = TRUE))

# Extract packages called using library()
lib_calls <- grep('^library', code, ignore.case=T, value=TRUE)
lib_calls <- gsub('^library[(]', '', lib_calls)
lib_calls <- gsub('[)]', '', lib_calls)
lib_calls <- gsub('^library$', '', lib_calls)
lib_calls <- gsub('["\']', '', lib_calls)  # Remove quotation marks

# Extract packages called using pkg_function_calls
pkg_function_calls <- regmatches(code, regexpr("\\b\\w+::", code))
pkg_function_calls <- gsub("::", "", unlist(pkg_function_calls))

# Combine and get unique package names
uniq_packages <- unique(c(lib_calls, pkg_function_calls))

# Remove any "empty" package names
uniq_packages <- uniq_packages[!uniq_packages == '']

# Organize the list alphabetically
uniq_packages <- uniq_packages[order(uniq_packages)]

cat('Required packages: \n')
cat(paste0(uniq_packages, collapse= ', '), fill=T)
cat('\n\n')

# Fetch the list of currently installed packages
installed_packages <- installed.packages()[, 'Package']

# Packages need to be installed
to_be_installed <- setdiff(uniq_packages, installed_packages)

# Check availability of packages for the current R version
available_packages <- rownames(
  available.packages(repos = 'https://cloud.r-project.org')
)
installing <- intersect(to_be_installed, available_packages)
not_available <- setdiff(to_be_installed, available_packages)


if (length(installing) > 0) {
  cat('Installing packages:\n')
  cat(paste0(installing, collapse= ', '), fill=T)
  cat('\n\n')
}

if (length(not_available) > 0) {
  cat('Packages not available on rcran:\n')
  cat(paste0(not_available, collapse= ', '), fill=T)
  cat('\n\n')
}

if (length(to_be_installed) == 0) {
  cat('All packages installed already!\n')
}

# Install the missing packages
if (length(installing) > 0) {
  install.packages(installing, type="source",
                   repos = 'https://cloud.r-project.org')
}

cat('\nDone!\n\n')
