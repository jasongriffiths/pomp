pompBuilder <- function (data, times, t0, name,
                         statenames, paramnames, tcovar, covar,
                         rmeasure, dmeasure, step.fn, step.fn.delta.t,
                         skeleton, skeleton.type = c("map","vectorfield"),
                         skelmap.delta.t = 1,
                         parameter.transform, parameter.inv.transform,
                         rprior, dprior,
                         globals, ..., link = TRUE, save = FALSE) {
  
  if (!is.data.frame(data)) stop(sQuote("data")," must be a data-frame")
  obsnames <- names(data)
  obsnames <- setdiff(obsnames,times)
  if (!missing(covar)) {
    if (!is.data.frame(covar)) stop(sQuote("covar")," must be a data-frame")
    covarnames <- colnames(covar)
    covarnames <- setdiff(covarnames,tcovar)
  } else {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
    covarnames <- character(0)
  }
  skeleton.type <- match.arg(skeleton.type)

  if (missing(statenames)) stop(sQuote("statenames")," must be supplied");
  if (missing(paramnames)) stop(sQuote("paramnames")," must be supplied");

  mpt <- missing(parameter.transform)
  mpit <- missing(parameter.inv.transform)
  if (xor(mpt,mpit))
    stop("if you supply one transformation function, you must supply its inverse")

  pompCBuilder(
               name=name,
               statenames=statenames,
               paramnames=paramnames,
               covarnames=covarnames,
               obsnames=obsnames,
               rmeasure=rmeasure,
               dmeasure=dmeasure,
               step.fn=step.fn,
               skeleton=skeleton,
               parameter.transform=parameter.transform,
               parameter.inv.transform=parameter.inv.transform,
               rprior=rprior,
               dprior=dprior,
               globals=globals,
               link=link,
               save=save
               ) -> name

  pomp(
       data=data,
       times=times,
       t0=t0,
       rprocess=euler.sim(
         step.fun=render(fnames$step.fn,name=name),
         delta.t=step.fn.delta.t,
         PACKAGE=name
         ),
       rmeasure=render(fnames$rmeasure,name=name),
       dmeasure=render(fnames$dmeasure,name=name),
       skeleton=render(fnames$skeleton,name=name),
       skeleton.type=skeleton.type,
       skelmap.delta.t=skelmap.delta.t,
       parameter.transform=render(fnames$parameter.transform,name=name),
       parameter.inv.transform=render(fnames$parameter.inv.transform,name=name),
       rprior=render(fnames$rprior,name=name),
       dprior=render(fnames$dprior,name=name),
       PACKAGE=name,
       statenames=statenames,
       paramnames=paramnames,
       tcovar=tcovar,
       covar=covar,
       ...
       )
}

pompLink <- function (name) {
  solib <- paste0(name,.Platform$dynlib.ext)
  dyn.load(solib)
}

pompUnlink <- function (name) {
  solib <- paste0(name,.Platform$dynlib.ext)
  dyn.unload(solib)
}

define <- list(
               var="#define {%variable%}\t({%ptr%}[{%ilist%}[{%index%}]])\n",
               var.alt="#define {%variable%}\t({%ptr%}[{%index%}])\n"
               )

undefine <- list(
                 var="#undef {%variable%}\n"
                 )

header <- list(
               file="/* pomp model file: {%name%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n",
               rmeasure="\nvoid {%name%}_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               dmeasure= "\nvoid {%name%}_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               step.fn="\nvoid {%name%}_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
               skeleton="\nvoid {%name%}_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               parameter.transform="\nvoid {%name%}_par_trans (double *__pt, double *__p, int *__parindex)\n{\n",
               parameter.inv.transform="\nvoid {%name%}_par_untrans (double *__pt, double *__p, int *__parindex)\n{\n",
               rprior="\nvoid {%name%}_rprior (double *__p, int *__parindex)\n{\n",
               dprior="\nvoid {%name%}_dprior (double *__lik, double *__p, int give_log, int *__parindex)\n{\n"
               )


fnames <- list(
               rmeasure="{%name%}_rmeasure",
               dmeasure= "{%name%}_dmeasure",
               step.fn="{%name%}_stepfn",
               skeleton="{%name%}_skelfn",
               parameter.transform="{%name%}_par_trans",
               parameter.inv.transform="{%name%}_par_untrans",
               rprior="{%name%}_rprior",
               dprior="{%name%}_dprior"
               )

decl <- list(
             periodic_bspline_basis_eval="\tvoid (*periodic_bspline_basis_eval)(double,double,int,int,double*);\nperiodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n",
             get_pomp_userdata_int="\tconst int * (*get_pomp_userdata_int)(const char *);\nget_pomp_userdata_int = (const int *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_int\");\n",
             get_pomp_userdata_double="\tconst double * (*get_pomp_userdata_double)(const char *);\nget_pomp_userdata_double = (const double *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_double\");\n",
             `get_pomp_userdata(\\b|[^_])`="\tconst SEXP (*get_pomp_userdata)(const char *);\nget_pomp_userdata = (const SEXP (*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata\");\n"
             )

footer <- list(
               rmeasure="\n}\n\n",
               dmeasure="\n}\n\n",
               step.fn="\n}\n\n",
               skeleton="\n}\n\n",
               parameter.transform="\n}\n\n",
               parameter.inv.transform="\n}\n\n",
               rprior="\n}\n\n",
               dprior="\n}\n\n",
               globals="\n"
               )

utility.fns <- list()

callable.decl <- function (code) {
  fns <- vapply(names(decl),grepl,logical(1),code,perl=TRUE)
  do.call(paste0,decl[fns])
}

missing.fun <- function (name) {
  paste0("  error(\"'",name,"' not defined\");")
}

pompCBuilder <- function (name, statenames, paramnames, covarnames, obsnames,
                          rmeasure, dmeasure, step.fn, skeleton,
                          parameter.transform, parameter.inv.transform,
                          rprior, dprior, globals, save = FALSE, link = TRUE)
{

  if (missing(name))
    name <- paste0("pomp",
                   paste(
                         format(
                                as.hexmode(ceiling(runif(n=2,max=2^24))),
                                upper.case=TRUE
                                ),
                         collapse=""
                         )
                   )

  has.trans <- !(missing(parameter.transform))

  if (missing(globals)) globals <- ""

  name <- cleanForC(name)
  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  stem <- if (save) name else file.path(tempdir(),name)
  modelfile <- paste0(stem,".c") 
  solib <- paste0(stem,.Platform$dynlib.ext)

  if (.Platform$OS.type=="unix") {
    pompheader <- "pomp.h"
  } else {
    pompheader <- system.file("include/pomp.h",package="pomp")
  }
  
  out <- file(description=modelfile,open="w")
  
  cat(file=out,render(header$file,name=name,pompheader=pompheader))

  for (f in utility.fns) {
    cat(file=out,f)
  }

  cat(file=out,globals,footer$globals)

  ## variable/parameter/observations definitions
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=v-1))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(define$var,variable=covarnames[v],ptr='__covars',ilist='__covindex',index=v-1))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(define$var,variable=obsnames[v],ptr='__y',ilist='__obsindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=paste0("D",statenames[v]),ptr='__f',ilist='__stateindex',index=v-1))
  }
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paste0("T",paramnames[v]),ptr='__pt',ilist='__parindex',index=v-1))
  }
  cat(file=out,render(define$var.alt,variable="lik",ptr='__lik',index=0))

  if (has.trans) {
    ## parameter transformation function
    cat(file=out,render(header$parameter.transform,name=name))
    cat(file=out,callable.decl(parameter.transform))
    cat(file=out,parameter.transform,footer$parameter.transform)
    ## inverse parameter transformation function
    cat(file=out,render(header$parameter.inv.transform,name=name))
    cat(file=out,callable.decl(parameter.inv.transform))
    cat(file=out,parameter.inv.transform,footer$parameter.inv.transform)
  }

  ## rmeasure function
  if (missing(rmeasure)) rmeasure <- missing.fun("rmeasure")
  cat(file=out,render(header$rmeasure,name=name),rmeasure,footer$rmeasure)

  ## dmeasure function
  if (missing(dmeasure)) dmeasure <- missing.fun("dmeasure")
  cat(file=out,render(header$dmeasure,name=name),dmeasure,footer$dmeasure)

  ## Euler step function
  if (missing(step.fn)) step.fn <- missing.fun("step.fn")
  cat(file=out,render(header$step.fn,name=name))
  cat(file=out,callable.decl(step.fn))
  cat(file=out,step.fn,footer$step.fn)

  ## skeleton function
  if (missing(skeleton)) skeleton <- missing.fun("skeleton")
  cat(file=out,render(header$skeleton,name=name))
  cat(file=out,callable.decl(skeleton))
  cat(file=out,skeleton,footer$skeleton)

  ## rprior function
  if (missing(rprior)) rprior <- missing.fun("rprior")
  cat(file=out,render(header$rprior,name=name),rprior,footer$rprior)

  ## dprior function
  if (missing(dprior)) dprior <- missing.fun("dprior")
  cat(file=out,render(header$dprior,name=name),dprior,footer$dprior)

  ## undefine variables
  for (v in seq_along(paramnames)) {
    cat(file=out,render(undefine$var,variable=paramnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=statenames[v]))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(undefine$var,variable=covarnames[v]))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(undefine$var,variable=obsnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=paste0("D",statenames[v])))
  }
  for (v in seq_along(paramnames)) {
    cat(file=out,render(undefine$var,variable=paste0("T",paramnames[v])))
  }
  close(out)

  cflags <- paste0("PKG_CFLAGS=\"",
                  Sys.getenv("PKG_CFLAGS"),
                  " -I",system.file("include",package="pomp"),"\"")

  rv <- system2(
                command=R.home("bin/R"),
                args=c("CMD","SHLIB","-o",solib,modelfile),
                env=cflags
                )
  if (rv!=0)
    stop("cannot compile shared-object library ",sQuote(solib))
  else
    cat("model codes written to",sQuote(modelfile),
        "\nlink to shared-object library",sQuote(solib),"\n")

  if (link) {
    if (save) {
      pompLink(name)
    } else {
      pompLink(file.path(tempdir(),name))
    }
  }

  invisible(name)
}

cleanForC <- function (text) {
  text <- as.character(text)
  text <- gsub("\\.","_",text)
  text
}

render <- function (template, ...) {
  vars=list(...)
  n <- sapply(vars,length)
  if (!all((n==max(n))|(n==1)))
    stop("incommensurate lengths of replacements")
  short <- which(n==1)
  n <- max(n)
  for (i in short) vars[[i]] <- rep(vars[[i]],n)
  
  retval <- vector(mode="list",length=n)
  for (i in seq_len(n)) {
    tpl <- template
    for (v in names(vars)) {
      src <- sprintf("\\{%%%s%%\\}",v)
      tgt <- vars[[v]][i]
      tpl <- gsub(src,tgt,tpl)
    }
    retval[[i]] <- tpl
  }
  do.call(paste0,retval)
}
