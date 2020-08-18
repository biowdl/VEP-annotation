version 1.0

# Copyright (c) 2018 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/vep.wdl" as vep
import "tasks/samtools.wdl" as samtools
import "tasks/picard.wdl" as picard
import "tasks/vt.wdl" as vt

workflow Annotation {
    input {
        File vcfFile
        File vcfIndex
        File referenceFasta
        File referenceFastaFai
        Array[File]? customFiles
        Array[File]? customFileIndices
        Array[String]? customFields
        String outputDir = "."
        String cacheDir
        String cacheVersion
        Map[String, String] dockerImages = {
        "vep": "quay.io/biocontainers/ensembl-vep:100.1--pl526hecc5488_0",
        "vt": "quay.io/biocontainers/vt:0.57721--hdf88d34_2"
        }
    }
    
    call vt.Normalize as normalize {
        input:
            dockerImage = dockerImages["vt"],
            inputVCF = vcfFile,
            inputVCFIndex = vcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            outputPath = outputDir + "/normalized-" + basename(vcfFile)
    }

    call samtools.TabixList as list {
        input: 
            inputFile = vcfFile,
            indexFile = vcfIndex
    }

    scatter(chromosome in list.chromosomes) {

        call vep.Customs as customString {
            input:
                customFiles = customFiles,
                customFileIndices = customFileIndices,
                customFields = customFields
        }

        call vep.Annotation as annotation {
            input:
                dockerImage = dockerImages["vep"],
                vcfFile = normalize.outputVcf,
                cacheDir = cacheDir,
                cacheVersion = cacheVersion,
                outputPath = outputDir + "/chromosomes/annotated-" + chromosome + ".vcf.gz",
                chromosome = chromosome,
                customs = customsString.out
        }
        
        call samtools.Tabix as tabix {
            input:
                inputFile = annotation.outputVcf,
                outputFilePath = annotation.outputVcf
        }
    }

    call picard.GatherVcfs as gatherVcf{
        input:
            inputVcfs = annotation.outputVcf,
            inputVcfIndexes = tabix.index,
            outputVcfPath = outputDir + "/annotated-" + basename(vcfFile)
    }

    output {
        File annotatedVcf = gatherVcf.outputVcf
    }
}
