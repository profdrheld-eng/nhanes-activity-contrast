#!/usr/bin/env Rscript
# A partial earlier export must be rejected before input loading or output writes.
arg <- sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root <- dirname(dirname(normalizePath(arg)))
cases <- list(plot_density.R=c('figure_1_two_regime_binned_source.csv',
  'figure_1_two_regime_summary.csv','Figure_1_two_regime_density.png'),
  plot_flow.R=c('participant_flow_sankey.png','participant_flow_checks.csv'))
for(script in names(cases))for(target in cases[[script]]) {
  work <- tempfile('output-guard-');dir.create(file.path(work,'displays'),recursive=TRUE)
  path <- file.path(work,'displays',target)
  marker <- charToRaw('Preserve this existing output')
  writeBin(marker,path)
  output <- suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),
    c(shQuote(file.path(root,'R',script)),shQuote(work)),stdout=TRUE,stderr=TRUE))
  stopifnot(!is.null(attr(output,'status')),attr(output,'status')!=0,
            any(grepl('Refusing to overwrite',output)),
            identical(readBin(path,'raw',n=length(marker)),marker))
  unlink(work,recursive=TRUE)
}
cat('PASS: figure exports protect each output of an incomplete previous run\n')
