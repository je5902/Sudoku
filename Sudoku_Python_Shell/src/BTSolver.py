import SudokuBoard
import Variable
import Domain
import Trail
import Constraint
import ConstraintNetwork
import time
import random

class BTSolver:

    # ==================================================================
    # Constructors
    # ==================================================================

    def __init__ ( self, gb, trail, val_sh, var_sh, cc ):
        self.network = ConstraintNetwork.ConstraintNetwork(gb)
        self.hassolution = False
        self.gameboard = gb
        self.trail = trail

        self.varHeuristics = var_sh
        self.valHeuristics = val_sh
        self.cChecks = cc

    # ==================================================================
    # Consistency Checks
    # ==================================================================

    # Basic consistency check, no propagation done
    def assignmentsCheck ( self ):
        for c in self.network.getConstraints():
            if not c.isConsistent():
                return False
        return True

    """
        Part 1 TODO: Implement the Forward Checking Heuristic

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        Note: remember to trail.push variables before you assign them
        Return: a tuple of a dictionary and a bool. The dictionary contains all MODIFIED variables, mapped to their MODIFIED domain.
                The bool is true if assignment is consistent, false otherwise.
    """
    def forwardChecking ( self ):
        
        # trailStack is empty, so call arcConsistency to remove values from domains based on given nums
        if len(self.trail.trailStack) == 0:
            self.arcConsistency()
            return ({}, True)
       
        modifiedVars = dict()
       
        # get most recently assigned variable by getting last variable in trailStack
        currVar = self.trail.trailStack[-1][0]
        
        # for all the neighbors of currVar
        for neighbor in self.network.getNeighborsOfVariable(currVar):
            
            # if the neighbor is changeable, not assigned, and has the value of currVar in its domain
            if neighbor.isChangeable and not neighbor.isAssigned() and neighbor.getDomain().contains(currVar.getAssignment()):
                
                # push neighbor and its current domain onto the trail
                self.trail.push(neighbor)
                
                # remove the value from its domain
                neighbor.removeValueFromDomain(currVar.getAssignment())
                
                # if the neighbor's new domain is 0, return false 
                if neighbor.domain.size() == 0:
                    return ({}, False)
                        
                modifiedVars.update({neighbor : neighbor.getValues()})
            
        return (modifiedVars, True)
    
    # =================================================================
	# Arc Consistency
	# =================================================================
    def arcConsistency( self ):
        assignedVars = []
        for c in self.network.constraints:
            for v in c.vars:
                if v.isAssigned():
                    assignedVars.append(v)
        while len(assignedVars) != 0:
            av = assignedVars.pop(0)
            for neighbor in self.network.getNeighborsOfVariable(av):
                if neighbor.isChangeable and not neighbor.isAssigned() and neighbor.getDomain().contains(av.getAssignment()):
                    neighbor.removeValueFromDomain(av.getAssignment())
                    if neighbor.domain.size() == 1:
                        neighbor.assignValue(neighbor.domain.values[0])
                        assignedVars.append(neighbor)

    
    """
        Part 2 TODO: Implement both of Norvig's Heuristics

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        (2) If a constraint has only one possible place for a value
            then put the value there.

        Note: remember to trail.push variables before you assign them
        Return: a pair of a dictionary and a bool. The dictionary contains all variables 
		        that were ASSIGNED during the whole NorvigCheck propagation, and mapped to the values that they were assigned.
                The bool is true if assignment is consistent, false otherwise.
    """
    def norvigCheck ( self ):
        # FIRST STRATEGY OF NORVIG IS FORWARD CHECKING
        modifiedDict, isConsistent = self.forwardChecking()
        
        if not isConsistent:
            return ({}, False) 
        
        assignedDict = {}
        N = self.gameboard.N
        
        for v in modifiedDict:
            for c in self.network.getConstraintsContainingVariable(v):
                counter = [0]*(N+1)
                for v2 in c.vars:
                    for val in v2.getValues():
                        counter[val] += 1
                for i in range(1, N+1):
                    if counter[i] == 1:
                        for v3 in c.vars:
                            if v3.isChangeable() and not v3.isAssigned() and v3.getDomain().contains(i):
                                self.trail.push(v3)
                                v3.assignValue(i)
                                
                                _, isConsistent = self.forwardChecking()
                                if not isConsistent:
                                    return ({}, False)
                                else:
                                    assignedDict.update({v3 : i})

        return (assignedDict, True)

    """
         Optional TODO: Implement your own advanced Constraint Propagation

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournCC ( self ):
        return False

    # ==================================================================
    # Variable Selectors
    # ==================================================================

    # Basic variable selector, returns first unassigned variable
    def getfirstUnassignedVariable ( self ):
        for v in self.network.variables:
            if not v.isAssigned():
                return v

        # Everything is assigned
        return None

    """
        Part 2 TODO: Implement the Minimum Remaining Value Heuristic

        Return: The unassigned variable with the smallest domain
    """
    def getMRV ( self ):
        # Get all unassigned variables
        unassignedVars = []
        for v in self.network.variables:
            if not v.isAssigned():
                unassignedVars.append(v)
                
        if len(unassignedVars) == 0:
            return None # all variables are assigned, so return None
        # Set smallestDomainVar and Val to that of first variable in all unassigned
        smallestDomainVar = unassignedVars[0]
        smallestDomainVal = smallestDomainVar.domain.size()
        
        # loop through all unassigned varaibles, updating what the smallest domain and var found are
        for uv in unassignedVars:
            if uv.domain.size() < smallestDomainVal:
                smallestDomainVar = uv
                smallestDomainVal = smallestDomainVar.domain.size()
        return smallestDomainVar
    
    def _getNumUnassignedNeighbors(self, v):
        '''Helper function that returns the nnumber if unassigned neighbors of a given variable'''
        num = 0
        for neighbor in self.network.getNeighborsOfVariable(v):
            if not neighbor.isAssigned():
                num += 1
        return num
    
    """
        Part 2 TODO: Implement the Minimum Remaining Value Heuristic
                       with Degree Heuristic as a Tie Breaker

        Return: The unassigned variable with the smallest domain and affecting the  most unassigned neighbors.
                If there are multiple variables that have the same smallest domain with the same number of unassigned neighbors, add them to the list of Variables.
                If there is only one variable, return the list of size 1 containing that variable.
    """
    def MRVwithTieBreaker ( self ):
        # Get all unassigned variables
        unassignedVars = []
        for v in self.network.variables:
            if not v.isAssigned():
                unassignedVars.append(v)
        if len(unassignedVars) == 0:
            return [None] # all variables are assigned, so return None
        # Set smallestDomainVar and Val to that of first variable in all unassigned
        varsToReturn = [unassignedVars[0]]
        smallestDomainVal = unassignedVars[0].domain.size()
        numUnassignedNeighbors = self._getNumUnassignedNeighbors(unassignedVars[0])
        
        # loop through all unassigned varaibles, updating what the smallest domain, then breaking ties by checking unassigned neighbors
        for uv in unassignedVars:
            if uv.domain.size() < smallestDomainVal:
                smallestDomainVal = uv.domain.size()
                numUnassignedNeighbors = self._getNumUnassignedNeighbors(uv)
                varsToReturn = [uv]
            elif uv.domain.size() == smallestDomainVal:
                if self._getNumUnassignedNeighbors(uv) > numUnassignedNeighbors:
                    numUnassignedNeighbors = self._getNumUnassignedNeighbors(uv)
                    varsToReturn = [uv]
                elif self._getNumUnassignedNeighbors(uv) == numUnassignedNeighbors:
                    # only append var if it has both smallest domain and most unassigned neighbors
                    varsToReturn.append(uv)
        return varsToReturn

    """
         Optional TODO: Implement your own advanced Variable Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVar ( self ):
        return None

    # ==================================================================
    # Value Selectors
    # ==================================================================

    # Default Value Ordering
    def getValuesInOrder ( self, v ):
        values = v.domain.values
        return sorted( values )

    """
        Part 2 TODO: Implement the Least Constraining Value Heuristic

        The Least constraining value is the one that will knock the least
        values out of it's neighbors domain.

        Return: A list of v's domain sorted by the LCV heuristic
                The LCV is first and the MCV is last
    """
    def getValuesLCVOrder ( self, v ):
        # Get values
        values = v.domain.values
        
        # Create a dictionary were each domain value is a key that has a value corresponding 
        # to the # of values it would knock out of v's neighbors
        valDict = {val: 0 for val in values} 
        
        # Loop through all of v's neighbors and values to fill dictionary
        for neighbor in self.network.getNeighborsOfVariable(v):
            for val in values:
                if val in neighbor.domain.values:
                    valDict[val] += 1
        # Return list ordered by the values of valdict
        return [val for (val,_) in sorted(valDict.items(),key=lambda x:x[1])]

    """
         Optional TODO: Implement your own advanced Value Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVal ( self, v ):
        return None

    # ==================================================================
    # Engine Functions
    # ==================================================================

    def solve ( self, time_left=600):
        if time_left <= 60:
            return -1

        start_time = time.time()
        if self.hassolution:
            return 0

        # Variable Selection
        v = self.selectNextVariable()

        # check if the assigment is complete
        if ( v == None ):
            # Success
            self.hassolution = True
            return 0

        # Attempt to assign a value
        for i in self.getNextValues( v ):

            # Store place in trail and push variable's state on trail
            self.trail.placeTrailMarker()
            self.trail.push( v )

            # Assign the value
            v.assignValue( i )

            # Propagate constraints, check consistency, recur
            if self.checkConsistency():
                elapsed_time = time.time() - start_time 
                new_start_time = time_left - elapsed_time
                if self.solve(time_left=new_start_time) == -1:
                    return -1
                
            # If this assignment succeeded, return
            if self.hassolution:
                return 0

            # Otherwise backtrack
            self.trail.undo()
        
        return 0

    def checkConsistency ( self ):
        if self.cChecks == "forwardChecking":
            return self.forwardChecking()[1]

        if self.cChecks == "norvigCheck":
            return self.norvigCheck()[1]

        if self.cChecks == "tournCC":
            return self.getTournCC()

        else:
            return self.assignmentsCheck()

    def selectNextVariable ( self ):
        if self.varHeuristics == "MinimumRemainingValue":
            return self.getMRV()

        if self.varHeuristics == "MRVwithTieBreaker":
            return self.MRVwithTieBreaker()[0]

        if self.varHeuristics == "tournVar":
            return self.getTournVar()

        else:
            return self.getfirstUnassignedVariable()

    def getNextValues ( self, v ):
        if self.valHeuristics == "LeastConstrainingValue":
            return self.getValuesLCVOrder( v )

        if self.valHeuristics == "tournVal":
            return self.getTournVal( v )

        else:
            return self.getValuesInOrder( v )

    def getSolution ( self ):
        return self.network.toSudokuBoard(self.gameboard.p, self.gameboard.q)
