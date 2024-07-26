library(XML)
library(seqinr)
library(optparse)

# Function to parse FASTA file
parse_fasta <- function(fasta_file) {
  sequences <- read.fasta(fasta_file,forceDNAtolower = F)
  sequences <- lapply(sequences, function(seq) {paste0(seq, collapse = "")})
  return(sequences)
}

# Function to parse dates file
parse_dates <- function(date_file) {
  dates <- read.table(date_file, header = TRUE, stringsAsFactors = FALSE)
  return(dates)
}

# Function to recursively replace strings in all nodes and attributes
replace_strings_in_node <- function(node, old_prefix, new_prefix) {
  # Replace in node text
  if (!is.null(xmlValue(node)) && grepl(old_prefix, xmlValue(node))) {
    xmlValue(node) <- gsub(old_prefix, new_prefix, xmlValue(node))
  }
  
  # Replace in node attributes
  attrs <- xmlAttrs(node)
  if (length(attrs) > 0) {
    for (attr_name in names(attrs)) {
      if (grepl(old_prefix, attrs[[attr_name]])) {
        xmlAttrs(node)[[attr_name]] <- gsub(old_prefix, new_prefix, attrs[[attr_name]])
      }
    }
  }
  
  # Recursively apply to child nodes
  children <- xmlChildren(node)
  if (length(children) > 0) {
    for (child in children) {
      replace_strings_in_node(child, old_prefix, new_prefix)
    }
  }
}


# Function to update XML with new sequences and dates
update_xml <- function(xml_file, sequences, dates, old_prefix, new_prefix ,output_file, invariant = T) {
  
  xml <- xmlTreeParse(xml_file, useInternalNodes = TRUE)
  root <- xmlRoot(xml)
  
  # Update sequences in XML
  sequence_nodes <- getNodeSet(root, "//sequence")
  sequence_node<- sequence_nodes[1][[1]]
  taxon <- xmlGetAttr(sequence_node, "taxon")
  seq <- xmlGetAttr(sequence_node, "value")
  seq_node<-as.character(saveXML(sequence_nodes[[1]]))
  parent_node<-getNodeSet(xml, "//data")[[1]]
  for (i in 1:length(sequences)){
    replace_node<-gsub(taxon,names(sequences)[i],seq_node)
    replace_node<-gsub(seq,sequences[i],replace_node,fixed = T)
    replace_node<-substr(replace_node, 2, nchar(replace_node) - 2)
    newnode<-newXMLNode(replace_node)
    addChildren(parent_node,newnode)
  }
  for (sequence_node in sequence_nodes) {
    removeNodes(sequence_nodes)
  }
    
  # Add invariant sites information
  if (invariant){ # Note: these proportions are specific to the nucleotide frequencies in H37Rv
    parent_node<-getNodeSet(xml, "//data")[[1]]
    xmlAttrs(parent_node)["id"]<-paste0(new_prefix,"original")
  
    As<-round((758552/4411532)*(4411532-length(unlist(strsplit(unlist(as.character(sequences[1])),"")))),0)
    Cs<-round((1449998/4411532)*(4411532-length(unlist(strsplit(unlist(as.character(sequences[1])),"")))),0)
    Gs<-round((1444614/4411532)*(4411532-length(unlist(strsplit(unlist(as.character(sequences[1])),"")))),0)
    Ts<-round((758368/4411532)*(4411532-length(unlist(strsplit(unlist(as.character(sequences[1])),"")))),0)
    newInvariant<-paste(As,Cs,Gs,Ts)
    
    new_node <- newXMLNode("data", 
                           attrs = c(id=new_prefix, 
                                     spec="FilteredAlignment", 
                                     filter="-", 
                                     data=paste0("@",new_prefix,"original"), 
                                     constantSiteWeights=newInvariant))
    
    addSibling(parent_node, new_node)
  }
  
  
  # Update taxa dates in XML
  trait_nodes <- getNodeSet(root, "//trait")
  trait_node<-trait_nodes[1][[1]]
  newDateString<-character()
  for (i in 1:length(sequences)){
    newDateString<-paste0(newDateString,names(sequences)[i],"=",dates[match(names(sequences)[i],dates[,1]),2],",")
  }
  newDateString<-substr(newDateString, 1, nchar(newDateString) - 1)
  xmlAttrs(trait_node)["value"] <- newDateString
  
  # Replace old prefix with base name of fasta file
  replace_strings_in_node(root, old_prefix, new_prefix)
  
  # Save updated XML to output file
  saveXML(root, file = output_file)
}

# Main function to handle command-line arguments and call update functions
main <- function() {
  option_list <- list(
    make_option(c("-x", "--xml"), type = "character", help = "Input BEAST2 XML file"),
    make_option(c("-f", "--fasta"), type = "character", help = "Input FASTA file with new sequences"),
    make_option(c("-d", "--dates"), type = "character", help = "Input text file with sample names and dates"),
    make_option(c("-p", "--old_prefix"), type = "character", help = "Old prefix for log, trees etc. files in XML"),
    make_option(c("-n", "--new_prefix"), type = "character", help = "New prefix for log, trees etc. files in XML"),
    make_option(c("-o", "--output"), type = "character", help = "Output XML file with updated sequences and dates"),
    make_option(c("-i", "--invariant"), type = "logical", help = "Add line for invariant site information per nucleotide to rescale output")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  
  sequences <- parse_fasta(opt$fasta)
  dates <- parse_dates(opt$dates)
  
  update_xml(opt$xml, sequences, dates, opt$old_prefix, opt$new_prefix, opt$output, opt$invariant)
}

if (!interactive()) {
  main()
}
