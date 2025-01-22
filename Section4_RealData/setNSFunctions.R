# Functions created by Duncan T. Lang to modify three functions in the search path
# to reduce processing time when fitting analysis models. 
# - For more information, see: https://github.com/duncantl/MJSRandPerm/blob/main/README.md

# Function used prior to fitting analysis models 
# - Format: 
  #   origFuns <- setNSFunctions()
setNSFunctions =
  function(funs = c(".get.outside.method", "stopifnot")) #, "[[.data.frame"))
  {
    oldFuns = list()
    
    if(".get.outside.method" %in% funs) {
      ns = getNamespace("emmeans")
      
      oldFuns$.get.outside.method = list(ns = ns, fun = get(".get.outside.method", ns))
      unlockBinding(".get.outside.method", ns)
      assign(".get.outside.method", function(generic, cls) NULL, ns)
      lockBinding(".get.outside.method", ns)
    }
    
    
    if("stopifnot" %in% funs) {
      ns = getNamespace("base")
      oldFuns$stopifnot = list(ns = ns, fun = stopifnot)
      unlockBinding("stopifnot", ns)
      assign("stopifnot", function (..., exprs, exprObject, local = TRUE) { list(...); if(!missing(exprs)) exprs; if(!missing(exprObject)) exprObject;  invisible()}, ns)
      lockBinding("stopifnot", ns)
    }
    
    if("[[.data.frame" %in% funs) {    
      ns = getNamespace("base")
      f = get("[[.data.frame", ns)
      oldFuns$"[[.data.frame" = list(ns = ns, fun = f)
      body(f) = body(f)[-3]
      
      unlockBinding("[[.data.frame", ns)
      assign("[[.data.frame", f, ns)
      lockBinding("[[.data.frame", ns)
    }
    
    
    invisible(oldFuns)
  }

# Function used to restore functions after fitting analysis models
# - Format:
  #   resetNSFunctions(origFuns)
resetNSFunctions =
  function(prev)
  {
    invisible( mapply(function(var, x)  {
      unlockBinding(var, x$ns)
      assign(var, x$fun, x$ns)
      lockBinding(var, x$ns)        
    },  names(prev), prev))
  }