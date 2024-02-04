# -*- coding: utf-8 -*-
import FreeCAD,Draft
import PySide
from PySide import QtGui ,QtCore
from PySide.QtGui import *
from PySide.QtCore import *

def sequential_Sport_check (portList):
    checkVector = []
    checkType = []
    modes = []

    i=0
    while (i < len(portList)):
        modes = portList[i].modes
        j=0;
        while (j < len(modes)):
            checkVector.append(modes[j].SportNumber)
            checkType.append(modes[j].type)
            j = j+1
        i = i+1

    # find the number of S-ports
    SportCount=0
    i=0
    while (i < len(checkVector)):
       if (checkVector[i] > SportCount):
           SportCount = checkVector[i]
       i = i+1

    # check for existence of each S-port
    Sport=1
    while (Sport <= SportCount):
       foundVoltage = False
       foundCurrent = False

       i=0
       while (i < len(checkVector)):
           if (checkVector[i] == Sport):
              if (checkType[i] == "voltage"):
                  if (foundVoltage):
                      App.Console.PrintError("ERROR: Duplicate voltage S-port definition found for S-port number " + str(Sport) + ".\n")
                      return True
                  foundVoltage = True
              if (checkType[i] == "current"):
                  if (foundCurrent):
                      App.Console.PrintError("ERROR: Duplicate current S-port definition found for S-port number " + str(Sport) + ".\n")
                      return True
                  foundCurrent = True
           i = i+1


       if (not foundVoltage and not foundCurrent):
           App.Console.PrintError("ERROR: Missing S-port number " + str(Sport) + ".\n")
           return True

       Sport = Sport+1

    return False

def get_commandSymbol (label,name):
    if name[0:4] == "Text":
        return label[0:2]
    return "null"

def get_commandName (label):
    splitText=label[2:].split('(')
    if (len(splitText) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"
    if (splitText[0] == ''):
        App.Console.PrintError("ERROR: \"" + label + "\" is missing a name.\n")
        return "null"
    return splitText[0]

def has_oneOption (label):
    splitText1=label.split('(')
    if (len(splitText1) != 2):
        return False

    splitText2=splitText1[1].split(')')
    if (len(splitText2) != 2):
        return False

    splitText3=splitText2[0].split(',')
    if (len(splitText3) != 1):
        return False

    return True

# 1 option - (option)
def get_commandOption1 (label):
    splitText1=label.split('(')
    if (len(splitText1) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText2=splitText1[1].split(')') 
    if (len(splitText2) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText3=splitText2[0].split(',')
    if (len(splitText3) != 1):
        print ("splitText2[0]=" + splitText2[0] + "\n")
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    return splitText3[0]

# 2 options - (option1,option2)
# select with index=1 or index=2
def get_commandOption2 (label,index):
    splitText1=label.split('(')
    if (len(splitText1) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText2=splitText1[1].split(')')  
    if (len(splitText2) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText3=splitText2[0].split(',')
    if (len(splitText3) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    return splitText3[index-1]

def get_paths (label):
    splitText1=label[2:].split('{')
    if (len(splitText1) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    splitText2=splitText1[1].split('}')
    if (len(splitText2) != 2):
        App.Console.PrintError("ERROR: \"" + label + "\" is incorrectly formatted.\n")
        return "null"

    paths=splitText2[0].split(',')
    if (len(paths) == 0):
        App.Console.PrintError("ERROR: \"" + label + "\" is missing paths.\n")
        return "null"

    return paths


class Path:
    def __init__(self, name, closed):
        self.name = name 
        self.closed = closed
        self.x = []
        self.y = []
        self.z = []

    def add_point(self,x,y,z):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

    def print(self,file):
        file.write("Path\n")
        file.write("   name=" + self.name + "\n")

        i=0
        while (i < len(self.x)):
            file.write("   point=(" + str(self.x[i]) + "," + str(self.y[i]) + "," + str(self.z[i]) + ")" + "\n")
            i=i+1

        file.write("   closed=" + self.closed + "\n")
        file.write("EndPath\n")
        file.write("\n")

class Boundary:
    def __init__(self, name, type, material, wave_impedance):
        self.name = name
        self.type = type
        self.material = material
        self.wave_impedance = wave_impedance
        self.paths = []
        self.directions = []

    def add_path(self,direction,path):
        self.directions.append(direction);
        self.paths.append(path)

    def print(self,file):
        file.write("Boundary\n")
        file.write("   name=" + self.name + "\n")
        file.write("   type=" + self.type + "\n")
        if self.type == "surface_impedance":
            file.write("   material=" + self.material + "\n")
        if self.type == "radiation":
            if self.wave_impedance != "":
                file.write("   wave_impedance=" + self.wave_impedance + "\n")

        i=0
        while (i < len(self.paths)):
            if i == 0:
                file.write("   path=" + self.directions[0] + self.paths[0] + "\n")
            else:
                file.write("   path" + self.directions[i] + "=" + self.paths[i] + "\n")
            i=i+1

        file.write("EndBoundary\n")
        file.write("\n")

    def check(self,pathList):
        fail = False

        for path in self.paths:
            found = False
            for testpath in pathList:
                if path == testpath.name:
                    found = True
            if not found:
                App.Console.PrintError("ERROR: Boundary " + self.name + " calls for the undefined path \"" + path + "\".\n")
                fail = True

        if self.type != "surface_impedance" and self.type != "perfect_electric_conductor" and self.type != "perfect_magnetic_conductor" and self.type != "radiation":
            App.Console.PrintError("ERROR: Boundary " + self.name + " calls for an invalid type [must be \"SI\", \"PEC\", \"PMC\", or \"radiation\"].\n")
            fail = True

        return fail


class Mode:
    def __init__(self, portName, SportNumber, type):
        self.portName = portName
        self.SportNumber = SportNumber
        self.type = type
        self.paths = []
        self.directions = []
        self.category = "modal"
        self.used = False

    def add_path(self,direction,path):
        self.directions.append(direction);
        self.paths.append(path)

    def set_category(self,category):
        self.category=category

    def print(self,file):
        if self.category == "modal":
            file.write("   Mode\n")
        if self.category == "line":
            file.write("   Line\n")
        file.write("      Sport=" + str(self.SportNumber) + "\n")
        file.write("      type=" + self.type + "\n")

        i=0
        while (i < len(self.paths)):
            if i == 0:
                file.write("      path=" + self.directions[0] + self.paths[0] + "\n")
            else:
                file.write("      path" + self.directions[i] + "=" + self.paths[i] + "\n")
            i=i+1

        if self.category == "modal":
            file.write("   EndMode\n")
        if self.category == "line":
            file.write("   EndLine\n")

    def dump(self):
        if self.category == "modal":
            print("Mode")
        if self.category == "line":
            print("Line")
        print("   port=" + str(self.portName))
        print("   Sport=" + str(self.SportNumber))
        print("   type=" + self.type)

        i=0
        while (i < len(self.paths)):
            if i == 0:
                print("   path=" + self.directions[0] + self.paths[0])
            else:
                print("   path" + self.directions[i] + "=" + self.paths[i])
            i=i+1

        if self.category == "modal":
            print("EndMode")
        if self.category == "line":
            print("EndLine")

    def check(self,pathList):
        fail = False

        if self.SportNumber < 1:
            App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" S-port "  + str(self.SportNumber) + " uses an invalid port number [must be > 0].\n")
            fail = True

        for path in self.paths:
            found = False
            for testpath in pathList:
                if path == testpath.name:
                    found = True
            if not found:
                App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" S-port " + str(self.SportNumber) + " calls for the undefined path \"" + path + "\".\n")
                fail = True

        if self.type != "voltage" and self.type != "current":
            App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" S-port " + str(self.SportNumber) + " calls for an invalid type [must be \"voltage\" or \"current\"].\n")
            fail = True

        return fail

class Port:
    def __init__(self, portName, impedance_definition):
        self.portName = portName
        self.impedance_definition = impedance_definition
        self.impedance_calculation = "modal"
        self.paths = []
        self.directions = []
        self.modes = []

    def set_impedance_calculation(self):
        for mode in self.modes:
            if mode.category == "line":
                self.impedance_calculation="line"

    def add_path(self,direction,path):
        self.directions.append(direction);
        self.paths.append(path)

    def add_mode(self,mode):
        self.modes.append(mode)

    def print(self,file):
        file.write("Port\n")
        file.write("   name=" + str(self.portName) + "\n")

        i=0
        while (i < len(self.paths)):
            if i == 0:
                file.write("   path=" + self.directions[0] + self.paths[0] + "\n")
            else:
                file.write("   path" + self.directions[i] + "=" + self.paths[i] + "\n")
            i=i+1

        file.write("   impedance_definition=" + self.impedance_definition + "\n")
        file.write("   impedance_calculation=" + self.impedance_calculation + "\n")

        for mode in self.modes:
            mode.print(file)

        file.write("EndPort\n")
        file.write("\n")

    def check(self,pathList):
        fail = False

        for path in self.paths:
            found = False
            for testpath in pathList:
                if path == testpath.name:
                    found = True
            if not found:
                App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" calls for the undefined path \"" + path + "\".\n")
                fail = True

        if self.impedance_definition != "VI" and self.impedance_definition != "PV" and self.impedance_definition != "PI":
            App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" calls for an invalid impedance definition [must be \"VI\", \"PV\", or \"PI\"].\n")
            fail = True

        if self.impedance_calculation != "modal" and self.impedance_calculation != "line":
            App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" calls for an invalid impedance calculation [must be \"modal\" or \"line\"].\n")
            fail = True

        for mode in self.modes: 
            if mode.check(pathList):
                fail = True

        i=0
        while (i < len(self.modes)-1):
            j=i+1
            while (j < len(self.modes)):
                if (self.modes[i].SportNumber == self.modes[j].SportNumber):
                    if (self.modes[i].type == self.modes[j].type):
                        App.Console.PrintError("ERROR: Port \"" + str(self.portName) + "\" specifies modes and/or lines with duplicate S-port numbers " + str(self.modes[i].SportNumber) + "\n")
                        fail = True
                j = j+1
            i = i+1

        return fail

pathList = []
boundaryList = []
modeList = []
portList = []


doc = FreeCAD.ActiveDocument
objs = FreeCAD.ActiveDocument.Objects

file = doc.Name+"_ports.txt"
path = "./" + file

try:
    SaveName = QFileDialog.getSaveFileName(None,QString.fromLocal8Bit("Save the ports file"),path,"*.txt") # PyQt4
except Exception:
    SaveName, Filter = PySide.QtGui.QFileDialog.getSaveFileName(None, "Save the ports file", path,"*.txt") # PySide

if SaveName == "":
    App.Console.PrintMessage("Process aborted.\n")
else:
    try:
        fail = False

        file = open(SaveName, 'w')
        try:

            for obj in objs:

                recognizedCommand = False

                #if (not obj.ViewObject.Visibility):
                #    continue

                name = obj.Name
                label = obj.Label

                commandSymbol=get_commandSymbol(label,name)

                if (label[0:2] == "_P"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    nameText=label[2:]
                    if (nameText == ''):
                       App.Console.PrintError("ERROR: \"" + label +"\" does not include a name.\n")
                       fail = True
                       continue

                    if (name[0:4] == "Line"):
                       newPath=Path(nameText,'false')
                       newPath.add_point(obj.Start.x, obj.Start.y, obj.Start.z)
                       newPath.add_point(obj.End.x, obj.End.y, obj.End.z)
                       pathList.append(newPath)

                    if (name[0:9] == "Rectangle"):
                       newPath=Path(nameText,'true')
 
                       x1 = obj.Placement.Base.x
                       y1 = obj.Placement.Base.y
                       z1 = obj.Placement.Base.z
                       newPath.add_point(x1,y1,z1)

                       rotationMatrix=obj.Placement.Rotation.toMatrix();

                       x2 = x1 + rotationMatrix.A11*obj.Length.Value
                       y2 = y1 + rotationMatrix.A21*obj.Length.Value
                       z2 = z1 + rotationMatrix.A31*obj.Length.Value
                       newPath.add_point(x2,y2,z2)

                       x3 = x1 + rotationMatrix.A11*obj.Length.Value + rotationMatrix.A12*obj.Height.Value
                       y3 = y1 + rotationMatrix.A21*obj.Length.Value + rotationMatrix.A22*obj.Height.Value
                       z3 = z1 + rotationMatrix.A31*obj.Length.Value + rotationMatrix.A32*obj.Height.Value
                       newPath.add_point(x3,y3,z3)

                       x4 = x1 + rotationMatrix.A12*obj.Height.Value
                       y4 = y1 + rotationMatrix.A22*obj.Height.Value
                       z4 = z1 + rotationMatrix.A32*obj.Height.Value
                       newPath.add_point(x4,y4,z4)

                       pathList.append(newPath)

                    if (name[0:4] == "Wire"):
                       closed="false"
                       if (obj.Closed):
                           closed="true"

                       x1 = obj.Placement.Base.x
                       y1 = obj.Placement.Base.y
                       z1 = obj.Placement.Base.z
                       rotationMatrix=obj.Placement.Rotation.toMatrix();

                       newPath=Path(nameText,closed)                       
                       for point in obj.Points:
                           x2 = x1 + rotationMatrix.A11*point.x + rotationMatrix.A12*point.y
                           y2 = y1 + rotationMatrix.A21*point.x + rotationMatrix.A22*point.y
                           z2 = z1 + rotationMatrix.A31*point.x + rotationMatrix.A32*point.y
                           newPath.add_point(x2,y2,z2)
                       pathList.append(newPath)

                if (commandSymbol == "_S"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    portName=get_commandName(label)
                    if (portName == "null"):
                        fail = True
                        continue

                    impedance=get_commandOption1(label)
                    if (impedance == "null"):
                        fail = True
                        continue

                    paths=get_paths(label)
                    if (paths == "null"):
                        fail = True
                        continue

                    newPort=Port(portName,impedance)

                    i=0
                    while (i < len(paths)):

                        direction="+"
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            direction="-"
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        newPort.add_path(direction,pathName)

                        i=i+1

                    portList.append(newPort)

                if (commandSymbol == "_B"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    portName=get_commandName(label)
                    if (portName == "null"):
                        fail = True
                        continue

                    type = ""
                    material = ""
                    wave_impedance = ""
                    if has_oneOption(label):
                        type=get_commandOption1(label)
                        if (type == "SI"):
                            App.Console.PrintError("ERROR: \"" + label + "\" must include the material parameter.\n")
                            fail = True
                        if (type == "null"):
                            fail = True
                            continue
                    else:
                        type=get_commandOption2(label,1)
                        if (type == "PEC" or type == "PMC"):
                            App.Console.PrintError("ERROR: \"" + label + "\" must not include a second parameter.\n")
                            fail = True
                        if (type == "null"):
                            fail = True
                            continue

                        material=get_commandOption2(label,2)
                        if (material == "null"):
                            fail = True
                            continue

                        wave_impedance=get_commandOption2(label,2)
                        if (material == "null"):
                            fail = True
                            continue

                    fullTypeName=""
                    if (type == "SI"): 
                        fullTypeName="surface_impedance"
                    if (type == "PEC"):
                        fullTypeName="perfect_electric_conductor"
                    if (type == "PMC"):
                        fullTypeName="perfect_magnetic_conductor"
                    if (type == "radiation"):
                        fullTypeName="radiation"

                    paths=get_paths(label)
                    if (paths == "null"):
                        fail = True
                        continue

                    newBoundary=Boundary(portName,fullTypeName,material,wave_impedance)

                    i=0
                    while (i < len(paths)):

                        direction="+"
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            direction="-"
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        newBoundary.add_path(direction,pathName)
                        
                        i=i+1

                    boundaryList.append(newBoundary)

                if (commandSymbol == "_M" or commandSymbol == "_L"):
                    print ("processing:" + label)
                    recognizedCommand = True

                    portName=get_commandName(label)
                    if (portName == "null"):
                        fail = True
                        continue

                    Sport=get_commandOption2(label,1)
                    if (Sport == "null"):
                       fail = True
                       continue

                    try:
                        SportNumber=int(Sport)
                    except Exception:
                        App.Console.PrintError("ERROR: \"" + label + "\" does not use an integer S_port number.\n")
                        fail = True
                        continue

                    type=get_commandOption2(label,2)
                    if (type == "null"):
                        fail = True
                        continue

                    paths=get_paths(label)
                    if (paths == "null"):
                        fail = True
                        continue

                    newMode=Mode(portName,SportNumber,type)
                    if (commandSymbol == "_M"):
                        newMode.set_category("modal")
                    if (commandSymbol == "_L"):
                        newMode.set_category("line")

                    i=0
                    while (i < len(paths)):

                        direction="+"
                        pathName=paths[i]

                        pathText=paths[i].split('-')
                        if (len(pathText) == 2):
                            direction="-"
                            pathName=pathText[1]

                        pathText=paths[i].split('+')
                        if (len(pathText) == 2):
                            pathName=pathText[1]

                        newMode.add_path(direction,pathName)

                        i=i+1

                    modeList.append(newMode)

                if not recognizedCommand and label[0:1] == '_':
                    App.Console.PrintError("ERROR: \"" + label + "\" calls for the invalid command \"" + commandSymbol + "\".\n")
                    fail = True

            # attach modes to ports
            for port in portList:
                for mode in modeList:
                    if port.portName == mode.portName:
                        port.add_mode(mode)
                        mode.used = True
                port.set_impedance_calculation()

            # warnings
            for mode in modeList:
               if (not mode.used):
                   App.Console.PrintError("Warning: Mode/Line on Port \"" + str(mode.portName) + "\" S-port " + str(mode.SportNumber) + " is not used within a port.\n")

            # error checks

            for boundary in boundaryList:
                if boundary.check(pathList):
                    fail=True

            for port in portList:
                if port.check(pathList):
                    fail=True

            i=0
            while (i < len(portList)-1):
                j=i+1
                while (j < len(portList)):
                    if (portList[i].portName == portList[j].portName):
                        App.Console.PrintError("ERROR: Port name \"" + str(portList[i].portName) + "\" is duplicated.\n")
                        fail = True
                    j = j+1
                i = i+1

            if sequential_Sport_check(portList):
                fail = True

            # write the file

            if not fail:

                # file type and version number
                file.write("#OpenParEMports 1.0\n\n")

                # File block
                file.write("File\n")
                file.write("   name=" + str(doc.FileName) + "\n")
                file.write("EndFile\n\n")

                # paths
                for path in pathList:
                    path.print(file)

                # boundaries
                for boundary in boundaryList:
                    boundary.print(file)

                # ports
                for port in portList:
                    port.print(file)


        except Exception:
            App.Console.PrintError("Undetermined fatal error while saving.\n")
        finally:
            file.close()
            if not fail:
                App.Console.PrintMessage("Saved file " + SaveName + ".\n")
    except Exception:
        App.Console.PrintError("Error Open file "+SaveName+"\n")


