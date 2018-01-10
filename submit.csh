#!/bin/csh

set prodId=`date +%F_%H-%M`

mkdir ./myOutput/${prodId}
mkdir ./jobs/log/${prodId}
mkdir ./jobs/err/${prodId}

star-submit-template -template submitToyMcZeroVtx.xml -entities productionId=${prodId}
