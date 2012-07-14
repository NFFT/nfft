<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
   <xsl:output method="xml" indent="yes" />

   <xsl:param name="suitename" />
   
   <xsl:variable name="cunitCount" select="count(//CUNIT_RUN_TEST_RECORD)"/>
   <xsl:variable name="cunitFailureCount" select="count(//CUNIT_RUN_TEST_FAILURE)"/>  

   <xsl:template match="/">
      <testsuites>
         <xsl:attribute name="errors">0</xsl:attribute>
         <xsl:attribute name="failures">
          <xsl:value-of select="$cunitFailureCount"/>
         </xsl:attribute>
         <xsl:attribute name="tests">
          <xsl:value-of select="$cunitCount"/>
         </xsl:attribute>
         <xsl:attribute name="name">
            <xsl:value-of select="$suitename" />
         </xsl:attribute>
         <xsl:apply-templates />
      </testsuites>
   </xsl:template>
   
  <xsl:template match="/CUNIT_TEST_RUN_REPORT/CUNIT_RESULT_LISTING">
     <xsl:for-each select="CUNIT_RUN_SUITE/CUNIT_RUN_SUITE_SUCCESS">
      <xsl:variable name="localCunitFailureCount" select="count(CUNIT_RUN_TEST_RECORD/CUNIT_RUN_TEST_FAILURE)"/>  
      <xsl:variable name="localCunitCount" select="count(CUNIT_RUN_TEST_RECORD)"/>
      <xsl:variable name="sn" select="normalize-space(SUITE_NAME/text())"/>
      <testsuite>
         <xsl:attribute name="errors">0</xsl:attribute>
         <xsl:attribute name="failures">
          <xsl:value-of select="$localCunitFailureCount"/>
         </xsl:attribute>
         <xsl:attribute name="tests">
          <xsl:value-of select="$localCunitCount"/>
         </xsl:attribute>
         <xsl:attribute name="name">
          <xsl:value-of select="$sn"/>
         </xsl:attribute>
         <xsl:apply-templates select="CUNIT_RUN_TEST_RECORD" />
      </testsuite>
     </xsl:for-each>
  </xsl:template>
      
  <xsl:template match="CUNIT_RUN_TEST_RECORD">
    <xsl:apply-templates select="CUNIT_RUN_TEST_SUCCESS" />
    <xsl:apply-templates select="CUNIT_RUN_TEST_FAILURE" />
  </xsl:template>
    
  <xsl:template match="CUNIT_RUN_TEST_SUCCESS">
    <testcase>
       <!--xsl:attribute name="classname">
          <xsl:value-of select="substring-before(substring-after(TEST_NAME,'test_'),'_')" />
       </xsl:attribute-->
       <xsl:attribute name="name">
          <xsl:value-of select="normalize-space(TEST_NAME)" />
       </xsl:attribute>
       <xsl:attribute name="time">0</xsl:attribute>
    </testcase>
  </xsl:template>
  <xsl:template match="CUNIT_RUN_TEST_FAILURE">
    <testcase>
       <!--xsl:attribute name="classname">
          <xsl:value-of select="substring-before(substring-after(TEST_NAME,'test_'),'_')" />
       </xsl:attribute-->
       <xsl:attribute name="name">
          <xsl:value-of select="normalize-space(TEST_NAME)" />
       </xsl:attribute>
       <xsl:attribute name="time">0</xsl:attribute>
       <failure>
         <xsl:attribute name="message">
            <xsl:value-of select=" normalize-space(CONDITION)" />
         </xsl:attribute>
         <xsl:attribute name="type">Failure</xsl:attribute>
         <xsl:value-of select="normalize-space(CONDITION)" />
         File: <xsl:value-of select="normalize-space(FILE_NAME)" />
         Line: <xsl:value-of select="normalize-space(LINE_NUMBER)" />
       </failure>
    </testcase>
  </xsl:template>
 
  <xsl:template match="text()|@*" />
  
</xsl:stylesheet>