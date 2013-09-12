#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
wiki2texi.py

Converts a mediawiki wiki formatted document to texinfo

Usage: wiki2texi.py PATH TITLE FILENAME
 PATH:  path to the wiki formatted file. 
 TITLE: texinfo title
 FILENAME: texinfo filename

The following wiki tags are not handled, conversion may unsuccessful when they are encountered:
    ; desc: item
    <small>
    <strike>
    <center>
    <!--comments-->
In contrast to mediawiki < > </ > style tags are not handled if they cross line boundaries.
"""

import re
import sys
from optparse import OptionParser  

def reverse(data):
    """
    Generator to reverse data
    """
    for index in range(len(data)-1, -1, -1):
        yield data[index]


class Wikilist:
    "Represents a list in the mediawiki document"
    
    endstring = {"*": "@end itemize", "#": "@end enumerate"}
    beginstring = {("*",0): "@itemize @bullet", ("*",1): "@itemize @minus", ("#",0): "@enumerate", ("#",1): "@enumerate"}
    def __init__(self,type,level):
        self.type = type
        self.level = level
    def __eq__(self,other):
        if isinstance(other,Wikilist):
            return(self.type == other.type and self.level == other.level)
        return False
    def __ne__(self,other):
        return(not self.__eq__(other))
    def begin(self):
        return(self.beginstring[(self.type,self.level)])
    def end(self):
        return(self.endstring[self.type])

class Listmanager:
    "Manages Wikilists"
    
    lists = []
    def query(self,symbols, listitem):
        newlist = Wikilist(symbols[0],len(symbols)-1)
        if not newlist in self.lists:
            self.lists.append(newlist)
            return(newlist.begin() +"\n@item " + listitem)
        if self.lists[-1] == newlist:
            return("@item " + listitem)
        else:
            ret = self.closelists(newlist)
            self.lists.pop()
            return(ret + "@item " + listitem + "\n" + newlist.end())
            
    def closelists(self,upto=None):
        ret = ""
        for list in reverse(self.lists):
            if list == upto:
                break
            ret += list.end() + "\n" 
            self.lists.pop()
        return ret

    def isopen():
        return(self.lists != [])


def sectionsub(matchobj):
    """
    lookup a section heading
    """
    switchsec = {"==": "@chapter ", "===": "@section ", "====": "@subsection ", "=====": "@subsubsection "}
    #texinfo supports subsubsection at max. More is relentlessly stripped.
    return(switchsec[matchobj.group(1)] + matchobj.group(2).strip("="))

def emphsub(matchobj):
    """
    lookup an emphasize command
    """
    switchemph = {"''": "@emph{", "'''": "@strong{", "'''''": "@emph{@strong{"}
    return(switchemph[matchobj.group(1)] + matchobj.group(2) + "}")


def convert(text,filename,title):
    """
    convert text from mediawiki syntax to texinfo
    """
    lists = Listmanager()
    verbopen = False
    
    #Create info Header 
    ret = " ".join(["\input texinfo  @c -*-texinfo-*-\n@setfilename", filename, "\n@settitle",  title, "\n\n@contents\n"])
    
    for line in text:
        #verbatim handling
        if not verbopen:
            #begin of a verbatim enviroment
            if re.match(" +\S+",line) != None:
                #start verbatim enviroment
                line = "\n@verbatim\n" + line[1:]
                #close open lists
                line = lists.closelists() +line
                verbopen = True
        else:
            #verbatim enviroment
            if re.match(" +\S+",line) != None:
                line = line[1:]
                # <> style commands work in mediawiki verbatim mode oO
                line = re.sub("<br>","\n",line)
                line = re.sub("<tt>(.+?)</tt>","\\1",line)
            else:
                line = "@end verbatim\n\n" + line
                verbopen = False
        if not verbopen:
            line = re.sub("\{\{i18n.*\}\}","",line) #strip i18n links
    
            line = re.sub("\[\[Category:.*\]\]","",line) #strip categories
            
            #escape @
            line = re.sub("@(?!end verbatim)","@@",line)           
            line = re.sub("\{\{ *Box Note *\|(.+?)\}\}","@cartouche\n\\1\n@end cartouche",line) #{{Box Note| }}

            #escape { and }
            line = line.replace("{","@{")
            line = line.replace("}","@}")

            line = re.sub("(={2,5})(.*)\\1",sectionsub,line) #section headings
    
            # <> command support (incomplete)
            line = re.sub("<br>","@\n",line) #<br>
            line = re.sub("<tt>(.+?)</tt>","@code{\\1}",line) #<tt>
            line = re.sub("<i>(.+?)</i>","@i{\\1}",line) #<i>
            line = re.sub("<b>(.+?)</b>","@b{\\1}",line) #<b>

            #Links
            line = re.sub("\[\[(.+?\|)?(.+?)\]\]","\\2",line) #Wiki link
            line = re.sub("\[(http://.+?) (.+?)\]","@uref{\\1, \\2}",line) #weblink
            line = re.sub("\[(https://.+?) (.+?)\]","@uref{\\1, \\2}",line) #weblink
            line = re.sub("\[(http://.+?)\]","@uref{\\1}",line) #weblink with description
            line = re.sub("\[(https://.+?)\]","@uref{\\1}",line) #weblink with description
    
            line = re.sub("('{2,3}|'{5})(.+?)\\1",emphsub,line) #'' and friends
    
            line = re.sub("(?m)^-{4}\s*",80*"-"+"@*",line)  #---- 

            #make : a single paragraph, we indent anyway
            line = re.sub("^:(.+)","\n\\1\n",line) 
            # Lists
            if (not (line.startswith("#") or line.startswith("*"))) and lists.isopen:
                line = lists.closelists() + line
            line = re.sub("(?m)^(\*+|#+)(.*)",lambda match: lists.query(match.group(1),match.group(2)),line)
            
            
        ret += line
    #close open enviroments
    if verbopen:
        ret += "@end verbatim\n@bye"
    else:
        ret += lists.closelists() +"@bye"
    return(ret)

def main(wikitext,filename,title):
    """
    Reads wiki formatted text in the file specified by argv or from stdin
    and prints texinfo to stdout

    Keyword Arguments:
     wikitext -- A list, containing lines of Mediawiki formatted text
     filename -- Filename for the texinfo file
     title -- Texinfo title
    """
    return((convert(wikitext,filename,title)))


if __name__ == '__main__':
	
    if len(sys.argv) < 4:
	use = '''
	 %prog  PATH  FILENAME TITLE
	 PATH:  path to the wiki formatted file. If ommited, input is read from stdin.
	 FILENAME: texinfo filename
	 TITLE: texinfo title
	'''
	parser = OptionParser(usage = use)
	(options, args) = parser.parse_args()
	parser.error("wrong number of arguments")
#        wikitext = sys.stdin.readlines()
#        print(main(wikitext,sys.argv[1],sys.argv[2]))
    else:
        wikitext = open(sys.argv[1]).readlines()
        print(main(wikitext,sys.argv[2],sys.argv[3]))

# vim: set ai ts=4 sw=4 et:
