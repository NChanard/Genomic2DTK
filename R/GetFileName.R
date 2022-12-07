#' Get file name.
#'
#' GetFileName
#' @description Function as `basename()` with the option to not return the file extension.
#' @param path.pth <character>: The path to the file.
#' @param ext.bln <logical>: Whether the file extension should be returned with the file name. (Default FALSE)
#' @return A character string.
#' @examples
#' filePath.pth ="my/path/to/my/file.txt"
#' GetFileName(path.pth=filePath.pth, ext.bln=FALSE)
#' GetFileName(path.pth=filePath.pth, ext.bln=TRUE)

GetFileName <- function(
    path.pth = NULL, ext.bln = FALSE
) {
    ifelse(ext.bln,
        fileName.str <- basename(path.pth),
        fileName.str <- sub(pattern = "(.*)\\..*$",
        replacement = "\\1", basename(path.pth))
    )
    return(fileName.str)
}