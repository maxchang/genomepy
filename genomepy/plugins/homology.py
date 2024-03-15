import re
import urllib.parse
import urllib.request

from loguru import logger

from genomepy.plugins import Plugin

ENSEMBL_PERMANENT_SITE = "apr2020.archive.ensembl.org"


class EnsemblHomologyPlugin(Plugin):
    homology_xml_query = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "CSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >

        <Dataset name = "{ensembl_shortname}_gene_ensembl" interface = "default" >
            <Filter name = "with_hsapiens_homolog" excluded = "0"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "hsapiens_homolog_ensembl_gene" />
            <Attribute name = "hsapiens_homolog_orthology_confidence" />
            <Attribute name = "hsapiens_homolog_orthology_type" />
        </Dataset>
    </Query>"""
    # <Dataset name = "${ensembl_shortname}_gene_ensembl" interface = "default" >


    def after_genome_download(self, genome, threads=1, force=False):
        fname = self.get_properties(genome)["human_homology"]
        # extract species name by parsing genome.readme_file
        ensembl_shortname = None
        with open(genome.readme_file, 'r') as readme:
            genome_filename = re.compile("(?<=original filename: )\\w+")
            for line in readme:
                species_match = genome_filename.search(line)
                # convert species name to ensembl short name
                if species_match:
                    tax_components = species_match[0].lower().split('_')
                    ensembl_shortname = ''.join([x[0] for x in tax_components[:-1]]) + tax_components[-1]
                    logger.info(f'For {genome.name}, using {ensembl_shortname} as Ensembl species short name.')
                    break
                # check provider
                elif line.startswith('provider:'):
                    if not 'Ensembl' in line:
                        logger.warning(f"{genome.name} was provdied by a source other than Ensembl.")
                        return

        if not ensembl_shortname:
            logger.warning(f"Could not parse species name for {genome.name}")
            return

        # insert correct organism
        query = urllib.parse.quote(self.homology_xml_query.format(ensembl_shortname=ensembl_shortname))
        url = f'http://{ENSEMBL_PERMANENT_SITE}/biomart/martservice?query={query}'
        try:
            urllib.request.urlretrieve(url, filename=fname)
        except Exception as e:
            logger.error(str(e))
            logger.error(f"Could not download homology file from {url}")


    def get_properties(self, genome):
        props = {"human_homology": re.sub(".fa(.gz)?$", ".human_homology.csv", genome.filename)}
        return props
