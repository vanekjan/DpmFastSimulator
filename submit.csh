#!/bin/csh

set prodId=`date +%F_%H-%M`

mkdir ./myOutput/${prodId}

star-submit-template -template submitPicoNpeAnaMaker.xml -entities productionId=${prodId}
