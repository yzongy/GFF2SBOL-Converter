
const fs = require('mz/fs')

const lines = (fs.readFileSync('yeast_chr11_3_34.gff') + '').split('\n')


const SBOLDocument = require('sboljs')

const uriPrefix = 'http://ncl.ac.uk/syntheticyeast/'

const annoPrefix = 'http://wiki.synbiohub.org/wiki/Terms/GFF3#'

const moment = require('moment')



const creationRegex = new RegExp(/# (.*) created from (.*) (.*) by (.*) \((.*)\)/g)
const subProvRegex = new RegExp(/# # (.*)/g)



var sbol = new SBOLDocument()

var SAs = {}

var rootCD = sbol.componentDefinition()
rootCD.addType('http://www.biopax.org/release/biopax-level3.owl#DnaRegion')
setIdentity(rootCD, uriPrefix, 'yeast_chr11_3_34')


var nAct = 0

var CDs = {
    [rootCD.uri]: rootCD
}


var inSequence = false


var nProv = countProv()
var curProv = null

var na = []

var skipFirstLine = true

lines.forEach((line) => {

    if(inSequence) {

        na.push(line)

    } else {

        if(line[0] === '#') {
            parseComment(line)
        } else if(line[0] === '>') {

            inSequence = true

        } else {
            parseRegularLine(line)
        }

    }


})

var seq = sbol.sequence()
setIdentity(seq, uriPrefix, 'yeast_chr11_3_34_seq')

seq.encoding = 'http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html'
seq.elements = na.join('').toLowerCase()

rootCD.addSequence(seq)


console.log(sbol.serializeXML({
    'xmlns:rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
    'xmlns:rdfs': 'http://www.w3.org/2000/01/rdf-schema#'
}))


function parseRegularLine(line) {

    /*if(skipFirstLine) {
        skipFirstLine = false
        return
    }*/

    const tokens = line.split('\t')

    const [
        seqid,
        source,
        type,
        start,
        end,
        score,
        strand,
        phase,
        attributes
    ] = tokens


    const attrMap = parseAttributes(attributes)


    var existingSA = SAs[attrMap.ID]

    if(existingSA) {

        addLocation(existingSA)


    } else {

        const sequenceAnnotation = sbol.sequenceAnnotation()
        rootCD.addSequenceAnnotation(sequenceAnnotation)

        setIdentity(sequenceAnnotation, rootCD.persistentIdentity + '/', attrMap.ID + '_anno')
        sequenceAnnotation.name = attrMap.Name


        SAs[attrMap.ID] = sequenceAnnotation

        addLocation(sequenceAnnotation)

        switch(type) {

            case 'chromosome':
            case 'CDS':
            case 'gene':


                const cd = createCD()

                const c = sbol.component()
                setIdentity(c, rootCD.persistentIdentity + '/', attrMap.ID + '_component')

                c.definition = cd

                sequenceAnnotation.component = c

                rootCD.addComponent(c)

                break

            default:

                break
        }

    }

    function addRolesAndMetadata(target) {

        const role = typeToRole(type)

        if(role)
            target.addRole(role)

        if(attrMap.Ontology_term) {

            attrMap.Ontology_term.split(',').forEach((term) => {
                if(term.indexOf('GO:') === 0) {
                    target.addRole('http://identifiers.org/go/' + term)
                }
            })
        }

        if(attrMap.Note) {

            target.description = decodeURIComponent(attrMap.Note)

        }

        target.addStringAnnotation(annoPrefix + 'source', source)
    }


    function parseAttributes(attributes) {

        const m = {}

        attributes.split(';').forEach((kv) => {

            const [ key, value ] = kv.split('=')

            m[key] = value
        })

        return m
    }

    function createCD() {

        const cd = sbol.componentDefinition()
        cd.addType('http://www.biopax.org/release/biopax-level3.owl#DnaRegion')

        setIdentity(cd, uriPrefix, attrMap.ID)
        cd.name = attrMap.Name

        addRolesAndMetadata(cd)

        return cd

    }

    function addLocation(sa) {

        const range = sbol.range()
        range.start = start
        range.end = end
        range.name = attrMap.Name

        setIdentity(range, sa.persistentIdentity + '/', 'range' + (sa._rangeN ? ++ sa._rangeN : (sa._rangeN = 1) && ''))

        if(strand === '+') {
            range.orientation = 'http://sbols.org/v2#inline'
        } else if(strand === '-') {
            range.orientation = 'http://sbols.org/v2#reverseComplement'
        } else if(strand === '.') {
        } else if(strand === '?') {
        }

        sa.addLocation(range)

    }

    function typeToRole(type) {

        switch(type) {

            case 'chromosome':
                return 'http://identifiers.org/so/SO:0000340'

            case 'CDS':
                return 'http://identifiers.org/so/SO:0000316'

            case 'telomere':
                return 'http://identifiers.org/so/SO:0000624'

            case 'deletion':
                return 'http://identifiers.org/so/SO:0000159'

            case 'site_specific_recombination_target_region':
                return 'http://identifiers.org/so/SO:0000342'

            case 'gene':
                return 'http://identifiers.org/so/SO:0000704'

            case 'stop_retained_variant':
                return 'http://identifiers.org/so/SO:0001567'

            case 'PCR_product':
                return 'http://identifiers.org/so/SO:0000006'

            case 'tag':
                return 'http://identifiers.org/so/SO:0000324'

            case 'restriction_enzyme_recognition_site':
                return 'http://identifiers.org/so/SO:0001687'

            case 'ARS':
                return 'http://identifiers.org/so/SO:0000436'

            case 'noncoding_exon':
                return 'http://identifiers.org/so/SO:0000198'

            default:
                break
        }

    }

}


function parseComment(line) {

    const r = creationRegex.exec(line)

    if(r) {

        const activity = sbol.provActivity()

        setIdentity(activity, uriPrefix, 'prov' + (nProv --))

        activity.name = r[5]
        activity.description = r[0].slice(2)
        activity.endedAtTime = moment(r[3], 'YYMMDD').toDate()

        const product = r[1]
        const productUri = uriPrefix + product + '/1'

        if(CDs[productUri]) {
            CDs[productUri].wasDerivedFrom = activity.uri
        } else {

            var cd = sbol.componentDefinition()
            cd.addType('http://www.biopax.org/release/biopax-level3.owl#DnaRegion')
            CDs[productUri] = cd

            setIdentity(cd, uriPrefix, product)

            cd.name = product
            cd.wasDerivedFrom = activity.uri
            cd.description = 'Automatically generated from GFF3 provenance; no further information available'
            cd.addUriAnnotation('http://www.w3.org/2000/01/rdf-schema#seeAlso', rootCD.uri)
        }

        const usage = sbol.provUsage()
        activity.addUsage(usage)

        setIdentity(activity, activity.persistentIdentity + '/', 'usage')

        usage.uri = activity.uri + '/usage'
        usage.entity = uriPrefix + r[2]
        usage.role = 'http://sbols.org/v2#source'

        curProv = activity

    } else {

        const r2 = subProvRegex.exec(line)

        if(r2) {

            curProv.description = curProv.description + '; ' + decodeURIComponent(r2[1])


        }

    }

}

function countProv() {

    var n = 0

    lines.forEach((line) => {

        if(creationRegex.exec(line))
            ++ n

    })

    return n
}

function setIdentity(target, prefix, displayId) {

    target.version = '1'
    target.displayId = cleanDisplayId(displayId)
    target.persistentIdentity = prefix + target.displayId
    target.uri = target.persistentIdentity + '/' + target.version

}

function cleanDisplayId(displayId) {

    return displayId.split('.').join('_')
                    .split('*').join('_')
                    .split(' ').join('_')
                    .split('(').join('_')
                    .split(')').join('_')
                    .split('-').join('_')
}



