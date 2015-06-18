#example of generation of search plan
#run this code, and then open with your web browser the file
fpath = "/tmp/plan.xml"

#more docs on XML
#http://127.0.0.1:19979/library/XML/html/newXMLDoc.html
#http://www.omegahat.org/RSXML/shortIntro.pdf


requiredPackages <- c("XML"
                      #,"doParallel"
                      )
packages <- requiredPackages[!(requiredPackages %in% rownames(installed.packages()))]
if(length(packages) > 0) {
  install.packages(packages)
}
invisible(lapply(requiredPackages, require, character.only=T))

#doc = newXMLDoc()
doc = newXMLDoc(dtd = "", namespaces=NULL, addFinalizer = TRUE, 
                name = "chiri search plan", node = NULL, isHTML = FALSE)
topNode = newXMLNode("top_node")
house1 = newXMLNode("house1", parent=topNode, attrs=c(address="Calle Grau 1", risk="High", lat=-10, long=30))
house2 = newXMLNode("house2", parent=topNode, attrs=c(address="Calle Grau 10", risk="High", lat=-10, long=30.1))

house1.1 = newXMLNode("house1_inspected__house4", parent=house1, attrs=c(address="Calle Franklin 42", risk="High", lat=-10, long=30.2))
house1.2 = newXMLNode("house1_inspected__house5", parent=house1, attrs=c(address="Calle Franklin 5", risk="High", lat=-10, long=30.21))
house1.3 = newXMLNode("house1_inspected__house6", parent=house1, attrs=c(address="Calle Franklin 6", risk="Medium", lat=-10, long=30.21))
house2.1 = newXMLNode("house2_inspected__house7", parent=house2, attrs=c(address="Calle Jefferson 2", risk="High", lat=-10, long=30.2))
#the node name indicates the branching structure, so that houseX_inspected is the condition, and now you go to houseY

inspectionPlanText <- saveXML(topNode)
cat(inspectionPlanText)

write(inspectionPlanText, file=fpath)
