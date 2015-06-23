#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.1

# history
# 25-Jun-2013 - original version


helptext="STATS DBSCAN VARIABLES=variable list
RDIST = number
METHOD = HYBRID* or RAW or DIST
MINPTS = number
/OPTIONS SCALE=YES or NO* SEEDS = YES* or NO
/OUTPUT PLOT=YES* or NO
/SAVE DSNAME=dataset name ID=variable name
RETAIN = YES or NO* WORKSPACE = file name
[/HELP]

This procedure performs density-based clustering.
The companion procedure, STATS DBPRED, can be used
to assign clusters to new data.

Example:
STATS DBSCAN VARIABLES = x y z RDIST = .5 MINPTS=5.

VARIABLES and RDIST are required.  All other parameters
are optional.

VARIABLES specifies the variables on which to cluster
the cases.  String variables are not permitted, and
categorical variables may be inappropriate for this
method.

Cases with any missing values are excluded listwise.

RDIST specifies the reachability value.  That is the
maximum Euclidean distance a point can be from another point to
be included in its cluster.  Some experimentation may
be required to find a satisfactory value.  Points more than
this distance from every other point are considered noise.
Too small a value will make every point noise while too large
a value will result in most points being in the same cluster.
Standardizing the data may make it easier to pick this value.

METHOD specifies the clustering method.  HYBRID is the default
and provides good speed with moderate memory usage.  RAW
uses less memory but is likely to be slow.  DIST clusters
based on a distance matrix.  For DIST, the input variables
must constitute a distance matrix.

MINPTS specifies the minimum number of points to constitute
a cluster.  The default value is 5.  It is recommended that
the value be at least the number of variables plus 1.  A
value of 2 is equivalent to a certain hierarchical
clustering algorithm.

SCALE = YES causes the variables to be standardized to zero
mean and unit variance.  This carries over to the prediction
process.

SEEDS = YES specifies that seed information is included in
the cluster object.  This is required if predicting new cases,
and statistics are reported for both seeds and border points
in the output when this value is YES.

PLOT=YES produces a plot of the points in the clusters.

The SAVE parameters control saving of the cluster information.
Specify DSNAME as a dataset name not already in use to save
the cluster assignments in a new dataset.  You can specify
an ID variable to identify the cases.

WORKSPACE specifies a file name for saving the R workspace.
RETAIN = YES causes the workspace to be retained in memory.
Either a workspace file or the retained workspace can be used
to predict new cases using STATS DBPRED.

STATS DBSCAN /HELP prints this help and does nothing else.
"

spssdbscan <- function(variables, rdist, method="hybrid",
    minpts=5, scale=FALSE, seeds=TRUE, dsname=NULL, id=NULL,
    retain=FALSE, workspace=NULL, plotit=TRUE, memorylimit=NULL) {

    setuplocalization("STATS_DBSCAN")

    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Density-Based Clustering")
    omsid="STATSDBSCAN"
    warns = Warn(procname=gtxt("Density-Based Clustering: Warnings"),omsid=omsid)

    tryCatch(library(fpc), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","fpc"),
            dostop=TRUE)
    }
    )
    if (!is.null(memorylimit)) {
        tryCatch(memory.limit(memorylimit), warning=function(e) {
            warns$warn(e, dostop=FALSE)
            } 
        )
    }
    
    # make sure active dataset has a name if creating any new ones
    # and names not in use
    alldsspecs = c(dsname)
    if (!is.null(alldsspecs)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name if creating new datasets"),
                dostop=TRUE)
        }
        if (length(intersect(alldsspecs, alldatasets) > 0)) {
            warns$warn(gtxt("One or more specified output dataset names are already in use"),
                dostop=TRUE)
        }
    }
    # The details list records all the settings
    details = list(variables=variables, rdist=rdist, method=method,
        minpts=minpts, scale=scale, dsname=dsname, seeds=seeds,id=id,
        retain=retain, workspace=workspace, plotit=plotit, 
        creationdate=as.character(Sys.time()))
    alldata = c(variables, id)
    nvars = length(variables)
    vardict = spssdictionary.GetDictionaryFromSPSS(alldata)
    if (max(apply(vardict["varType",1:nvars],1,as.integer)) > 0) {
        warns$warn(gtxt("String variables cannot be used in this procedure"),
            dostop=TRUE)
    }

    # dbscan does not support factors or missing data.
    # Any categorical variables are treated as numeric.
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE)
    allcasecount = nrow(dta)
    dta = dta[complete.cases(dta),]
    completecasecount = nrow(dta)
    if (!is.null(id)) {
        iddata= dta[nvars+1]
        dta = dta[-(nvars+1)]
    } else {
        iddata = NULL
    }
    
    # do the clustering
    res = tryCatch(dbscan(
        data=dta,
        eps=rdist,
        MinPts=minpts,
        scale=scale,
        method=method,
        seeds=seeds), error = function(e) {
            warns$warn(e$message, dostop=TRUE)
            }
        )
        

    # print results
    displayresults(warns, res, dta, details, allcasecount, completecasecount)
    warns$display(inproc=TRUE)
    createdataset(res, dta, iddata, details, vardict)

    # save and clean up workspace as needed
    # workspace can be used for prediction via STATS DBPRED

    if (retain || !is.null(workspace)) {
        assign("stats_dbscan_res", res, envir=.GlobalEnv)
        assign("stats_dbscan_details", details, envir=.GlobalEnv)
        assign("stats_dbscan_dta", dta, envir=.GlobalEnv)
        rm(res, details, dta)
        if (!is.null(workspace)) {
            save(stats_dbscan_res, stats_dbscan_details, stats_dbscan_dta, file=workspace)
        }
    }
    if (!retain) {
        res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
    }
}

displayresults = function(warns, res, dta, details, allcasecount, completecasecount) {
    # Display results tables and plot
    
    StartProcedure(gtxt("Density-Based Clustering"), "STATSDBSCAN")

    # summary results
    scaption = gtxt("Computations done by R package fpc function dbscan by Christian Hennig")
    lbls = c(
        gtxt("Variables"), 
        gtxt("Method"),
        gtxt("Min Reachability Distance"),
        gtxt("Data Scaled"),
        gtxt("Minimum Cluster Size"),
        gtxt("Output Dataset"), 
        gtxt("Workspace File"),
        gtxt("Workspace Retained"),
        gtxt("Creation Date"),
        gtxt("Number of Cases"),
        gtxt("Number of Valid Cases")
    )
    vals = c(
        paste(details[["variables"]], collapse=" "),
        details["method"],
        details["rdist"],
        ifelse(details["scale"], gtxt("Yes"), gtxt("No")),
        details["minpts"],
        ifelse(!is.null(details[["dsname"]]), details[["dsname"]], gtxt("--NA--")),
        ifelse(!is.null(details[["workspace"]]), details[["workspace"]], gtxt("--NA--")),
        ifelse(details[['retain']], gtxt("Yes"), gtxt("No")),
        details[["creationdate"]],
        allcasecount,
        completecasecount
    )

    # settings and result summary
    spsspivottable.Display(
        data.frame(cbind(vals), row.names=lbls), 
        title = gtxt("Settings and Results Summary"),
        collabels=c(gtxt("Summary")), templateName="DBSCANSUMMARY", outline=gtxt("Summary"),
        caption = scaption
    )

    # cluster distribution - just counts or counts for seed and border points separately
    if (is.null(res$isseed)) {
        clusterdist = data.frame(table(res$cluster))
        clabels = list(gtxt("Frequency"))
        clusterdist[1] = NULL
    } else {
        tt = table(res$cluster, res$isseed)
        if (ncol(tt) == 1) { # all points in the same category :-(
            if (dimnames(tt)[[2]] == "TRUE") {  #all points are seeds
                tt = cbind(rep(0, nrow(tt)), tt)  # isnert implied border counts
            } else {  # all points are borders
                tt = cbind(tt, rep(0, nrow(tt)))
            }
        }
        clusterdist = data.frame(seed=tt[,2], border=tt[,1])
        clusterdist[3] = apply(clusterdist, 1, sum)
        clabels = list(gtxt("Seed Frequency"), gtxt("Border Frequency"), gtxt("Total"))
    }
    clusterdist = rbind(clusterdist, colSums(clusterdist))
    row.names(clusterdist)[nrow(clusterdist)] = gtxt("Total")
    spsspivottable.Display(
        clusterdist, 
        title=gtxt("Cluster Membership Distribution"),
        rowdim=gtxt("Cluster Number"), hiderowdimtitle=FALSE,
        collabels=clabels,
        templateName="DBSCANDIST", 
        outline=gtxt("Distribution"),
        format = formatSpec.Count
    )

    if (details[["plotit"]] && nrow(clusterdist) > 1) {
        plot(res, as.matrix(dta), main=gtxt("Cluster Plot"))
    }
}

createdataset = function(res, dta, iddata, details, vardict) {
    # Create classification dataset if requested
    # Dataset name is known to be okay, and procedure state is ended

    # for existing data, the prediction is just
    # res$cluster
    # If eps is too small, all predictions will be 0 - noise
    
    if (!is.null(details[["dsname"]])) {
        pred = data.frame(res$cluster)
        if (!is.null(iddata)) {  # was an id variable provided
            pred = data.frame(iddata, pred)
            loc = match(details[["id"]], vardict["varName",])
            idtype = as.integer(vardict[["varType", loc]])
            idformat = vardict[["varFormat", loc]]
            idlabel = vardict[["varLabel", loc]]
        } else {
            pred = data.frame(id=1:nrow(pred), pred)
            idtype = 0
            idformat = "F10.0"
            idlabel = ""
        }

        dict = spssdictionary.CreateSPSSDictionary(
            c("ID", idlabel, idtype, idformat, "nominal"),
            c("PredCluster", gtxt("Predicted Cluster"), 0, "F8.0", "nominal")
        )
        spssdictionary.SetDictionaryToSPSS(details[["dsname"]], dict)
        spssdata.SetDataToSPSS(details[["dsname"]], pred)
        spssdictionary.EndDataStep()
    }
}
    
    
# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_DBSCAN"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_DBSCAN"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

    if (is.null(msg) || dostop) {
        lcl$display(inproc)  # display messages and end procedure state
        if (dostop) {
            stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
        }
    }
}

    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

    if (lcl$msgnum == 0) {   # nothing to display
        if (inproc) {
            spsspkg.EndProcedure()
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        } else {
            procok=TRUE
        }

        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
            rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

    for (i in 1:lcl$msgnum) {
        rowcategory = spss.CellText.String(as.character(i))
        BasePivotTable.SetCategories(table,rowdim,rowcategory)
        BasePivotTable.SetCellValue(table,rowcategory, 
            spss.CellText.String(lcl$msglist[[i]]))
    }
    spsspkg.EndProcedure()   # implies display
} else { # can't produce a table
    for (i in 1:lcl$msgnum) {
        print(lcl$msglist[[i]])
    }
}
}
}
return(lcl)
}
Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("VARIABLES", subc="",  ktype="existingvarlist", 
            var="variables", islist=TRUE),
        spsspkg.Template("RDIST", subc="",  ktype="float", var="rdist"),
        spsspkg.Template("METHOD", subc="",  ktype="str", 
            var="method", vallist=list("hybrid", "raw", "dist")),
        spsspkg.Template("MINPTS", subc="", ktype="int", var="minpts"),
        
        spsspkg.Template("SCALE", subc="OPTIONS",  ktype="bool", var="scale"),
        spsspkg.Template("SEEDS", subc="OPTIONS",  ktype="bool", var="seeds"),
        spsspkg.Template("MEMORYLIMIT", subc="OPTIONS", ktype="float",
            var="memorylimit", vallist=list(2047)),
        
        spsspkg.Template("DSNAME", subc="SAVE", ktype="varname", var="dsname"),
        spsspkg.Template("ID", subc="SAVE", ktype="existingvarlist", var="id"),
        spsspkg.Template("RETAIN", subc="SAVE", ktype="bool", var="retain"),
        spsspkg.Template("WORKSPACE", subc="SAVE", ktype="literal", var="workspace"),

        spsspkg.Template("PLOT", subc="OUTPUT", ktype="bool", var="plotit")
    ))        
if ("HELP" %in% attr(args,"names"))
    #writeLines(helptext)
    helper(cmdname)
else
    res <- spsspkg.processcmd(oobj,args,"spssdbscan")
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}