# loadModule.R
#
loadGitHub <- function( FuncTions ){
    library(RCurl)
    funcNum <- length(FuncTions)
    for( index in 1:funcNum){
        URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncTions[index])
        eval(parse(text = getURL(URL, ssl.verifypeer = FALSE)))
        cat(sprintf('Load %s  from GitHub.', URL)); cat('\n')
    }
    cat('Succeeded to load modules from GitHub.'); cat('\n')
    return()
}

loadLocal <- function( RPATH, FuncTions ){
    funcNum <- length(FuncTions)
    for( index in 1:funcNum){
        source( sprintf("%s/%s.R", RPATH, FuncTions[index]))
    }
    cat('Failed to access GitHub .... Load loal modules.'); cat('\n')
    return()
}
