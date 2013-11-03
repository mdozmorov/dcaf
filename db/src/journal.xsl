<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" />

<xsl:template match="MedlineCitation">
  <xsl:value-of select="translate(MedlineJournalInfo/NlmUniqueID, 'RA', '')"/>
  <xsl:text>&#x9;</xsl:text>
  <xsl:value-of select="MedlineJournalInfo/ISSNLinking"/>
  <xsl:text>&#x9;</xsl:text>
  <xsl:value-of select="MedlineJournalInfo/MedlineTA"/>
</xsl:template>

</xsl:stylesheet>
