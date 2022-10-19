<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet
	xmlns:p="http://www.predictprotein.org/predictprotein"
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform' 
	version='1.0'>
<xsl:output method="text"/>

<xsl:template match="/">
<xsl:for-each select="//p:feature">
<xsl:value-of select="/p:predictprotein/p:entry/p:accession"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="@type"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="p:location/p:begin/@position"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="p:location/p:end/@position"/>
<xsl:text>
</xsl:text>
</xsl:for-each>
<xsl:value-of select="//p:secondaryStructures/p:featureString"/>
<xsl:text>
</xsl:text>
</xsl:template>
</xsl:stylesheet>



