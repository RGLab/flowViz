context("ncdfFlowSet plots")
library(flowCore)
library(flowStats)
library(ncdfFlow)
data(GvHD)

morphGate <- norm2Filter("FSC-H", "SSC-H", filterId = "MorphologyGate",scale = 2)
fs <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
suppressMessages(ncfs <- ncdfFlowSet(fs))
samples <- sampleNames(ncfs)

test_that("xyplot ncdfFlowSet", {
  sn <- samples[1:2]
  nc <- ncfs[sn]
  ncObj <- xyplot(`SSC-H`~`FSC-H`, nc)
  
  expect_is(ncObj[["panel.args.common"]][["frames"]], "ncdfFlowSet")
  
  expect_equal(ncObj[["panel.args.common"]][["type"]], "xyplot")
  
  fsObj <- xyplot(`SSC-H`~`FSC-H`, fs[sn])
  
  expect_is(fsObj[["panel.args.common"]][["frames"]], "environment")
  
  expect_equal(sub("type", "plotType", deparse(ncObj[["call"]])), deparse(fsObj[["call"]]))
  
  ncObj[["panel.args.common"]][["frames"]] <- NULL
  fsObj[["panel.args.common"]][["frames"]] <- NULL
  ncObj[["panel.args.common"]][["type"]] <- NULL
  ncObj[["call"]] <- NULL
  fsObj[["call"]] <- NULL
  
  expect_equivalent(fsObj, ncObj)
})

test_that("densityplot ncdfFlowSet", {
  sn <- samples[1:2]
  nc <- ncfs[sn]
  ncObj <- densityplot(~`SSC-H`, nc)
  fsObj <- densityplot(~`SSC-H`, fs[sn])
  
  expect_is(ncObj[["panel.args.common"]][["frames"]], "ncdfFlowSet")
  expect_is(fsObj[["panel.args.common"]][["frames"]], "environment")
  
  expect_equal(deparse(ncObj[["call"]]), deparse(fsObj[["call"]]))
  
  ncObj[["panel.args.common"]][["frames"]] <- NULL
  fsObj[["panel.args.common"]][["frames"]] <- NULL
  ncObj[["call"]] <- NULL
  fsObj[["call"]] <- NULL
  
  expect_equivalent(fsObj, ncObj)
})