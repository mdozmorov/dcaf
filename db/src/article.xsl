<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" />

<xsl:template match="MedlineCitation">
    <xsl:value-of select="PMID" />
    <xsl:text>&#x9;</xsl:text>
    <xsl:value-of select="translate(MedlineJournalInfo/NlmUniqueID, 'AR', '')"/>
    <xsl:text>&#x9;</xsl:text>
    <xsl:choose>
        <xsl:when test="Article/Journal/JournalIssue/PubDate/Year">
            <xsl:value-of select="Article/Journal/JournalIssue/PubDate/Year" />
            <xsl:text>-</xsl:text>
            <xsl:choose>
                <xsl:when test="Article/Journal/JournalIssue/PubDate/Month">
                  <xsl:value-of select="Article/Journal/JournalIssue/PubDate/Month" />
                </xsl:when>
                <xsl:otherwise>01</xsl:otherwise>
            </xsl:choose>
            <xsl:text>-</xsl:text>
            <xsl:choose>
                <xsl:when test="Article/Journal/JournalIssue/PubDate/Day">
                  <xsl:value-of select="Article/Journal/JournalIssue/PubDate/Day" />
                </xsl:when>
                <xsl:otherwise>01</xsl:otherwise>
            </xsl:choose>
        </xsl:when>
    </xsl:choose>
    <xsl:text>&#x9;</xsl:text>
    <xsl:value-of select="translate(Article/ArticleTitle, '&#x9;&#xa;', ' ')" />
    <xsl:text>&#x9;</xsl:text>
    <xsl:value-of select="translate(Article/Abstract/AbstractText, '&#x9;&#xa;', ' ')" />
</xsl:template>

</xsl:stylesheet>
