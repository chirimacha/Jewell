"
This code will send an email from Amazon's AWS node

"

#TODO: install.packages("mailR")
library("mailR")

#Help on email
#https://aws.amazon.com/ec2/faqs/ 

from_addresses = c("ebillig@mail.med.upenn.edu")
to_addresses   = c("ebillig@mail.med.upenn.edu", "Sarah.Nutman@uphs.upenn.edu", "javierequintanilla@gmail.com")  #and so on
attachment_file_paths = c("/tmp/results_2015_07_23.csv")  #FIXME

send.mail(from=from_addresses, to=to_addresses, 
         subject = paste0("Results - Automated email", date()), 
         body = "Attached are the results from last night on Amazon", 
         attach.files = attachment_file_paths, 
         encoding = "iso-8859-1", html = FALSE, inline = FALSE, smtp = list(host.name="localhost"), authenticate = FALSE, send = TRUE, debug = FALSE)
