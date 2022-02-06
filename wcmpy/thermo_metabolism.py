# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 00:46:36 2022

@author: Andrew Freiburger
"""

class ThermoMet():
    def absorption(self,):
        #calculates the quantity of transferred electrons in the membrane
        membrane_reaction_quantity = 10000
        n_electrons_displaced = 0
        reaction_count = 0        
        for reaction in self.reactions:
            if reaction["compartment"] == "m":
                n_electrons_displaced += reaction["electrons per membrane reaction"] #!!! amendment 25
                reaction_count += 1
                  
        #calculates the ideal displacement from thermodynamic equilibrium
        average_electrons_displaced = n_electrons_displaced / reaction_count
        Q_keq_ratio = exp(-average_electrons_displaced*self.F*membrane_electrode_potential/R/optimum_incubation_temperature)
        
        #calculates the quantity of absorbed chemicals
        for potentially_absorbable_chemical in self.bulk_chemicals.iterrows():
            potentially_absorbed_chemicals["name"] = potentially_absorbable_chemical["name"]
            potentially_absorbed_chemicals["amount"] = self.bulk_chemicals.loc["name" == potentially_absorbable_chemical["name"], "concentration"] * volume_absorption_shell  
        
            #reactions with potentially absorbable chemicals are evaluated
            for index, row in self.reactions.iterrows():
                if row["Stoichiometry"] != None:
                    for chemical in row["Stoichiometry"]:
                        if chemical[self.index_name] == potentially_absorbable_chemical["name"]:
                            Q_numerator = 1
                            Q_denominator = 1
                            current_Q = 0
                            current_Q_keq = 0
                            current_keq_displacement = 0
                            
                            if chemical[self.index_compartment] == "c":     #!!! amendment 21 
                                if chemical[self.index_stoich_coef] > 0:  
                                    Q_numerator *= pow(cytoplasm_metabolites.at[chemical[self.index_name], "amount"], chemical[self.index_stoich_coef])                          
                                elif chemical[self.index_stoich_coef] < 0:    
                                    Q_denominator *= pow(cytoplasm_metabolites.at[chemical[self.index_name], "amount"], -1*chemical[self.index_stoich_coef])                            
    
                            elif chemical[self.index_compartment] == "m":
                                if chemical[self.index_stoich_coef] > 0:
                                    Q_numerator *= pow(self.cytoplasm_metabolites.at[chemical[self.index_name], "amount"], chemical[self.index_stoich_coef])
                                elif chemical[self.index_stoich_coef] < 0:
                                    Q_denominator *= pow(self.cytoplasm_metabolites.at[chemical[self.index_name], "amount"], -1*chemical[self.index_stoich_coef])
                                
                            if Q_denominator != 0:
                                current_Q = Q_numerator/Q_denominator
                                
                            if row["keq"] != "None":
                                if current_Q <= 1:
                                    current_Q_keq *= current_Q/float(row["keq"])  
                                    current_keq_displacement = Q_keq_ratio - current_Q_keq
            
            #absorbed chemicals is allocated into the cytoplasm
            if current_keq_displacement <= 1:
                A_absorbed = potentially_absorbed_chemicals[potentially_absorbable_chemical["name"], "amount"]*absorption_proportion*(current_keq_displacement)*potentially_absorbable_chemical["barrier"]
                bulk_chemicals.loc[potentially_absorbable_chemical["name"], "concentration"] -= A_absorbed   
                
                if potentially_absorbable_chemical["chemical_type"] == "antibiotic":
                    antibiotics.loc[potentially_absorbable_chemical["name"], "amount"] += A_absorbed
                    
                elif potentially_absorbable_chemical["chemical_type"] == "membrane_metabolite":
                    self.cytoplasm_metabolites.loc[potentially_absorbable_chemical["name"], "amount"] += A_absorbed
                    
                elif potentially_absorbable_chemical["chemical_type"] == "cytoplasm_metabolite":
                    cytoplasm_metabolites.loc[potentially_absorbable_chemical["name"], "amount"] += A_absorbed
                    
                    
    def possible_reactions_forward(self,reaction_stoichiometry):
        if str(reaction_stoichiometry) != "None":
            minimum = -1
            for metabolite in reaction_stoichiometry:   
                
                #the quantity of forward reaction progressions is determined
                if metabolite[self.index_stoich_coef] < 0:
                    if metabolite[self.index_compartment] == "c":
                        try:
                            temp = (self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"]/abs(metabolite[self.index_stoich_coef])).iloc[0]  
                            if temp < minimum:  
                                minimum = temp
                        except:
                            minimum = 0          
                            
                    elif metabolite[self.index_compartment] == "m":
                        try:
                            temp = (self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"]/abs(metabolite[self.index_stoich_coef])).iloc[0]
                            if temp < minimum:
                                minimum = temp
                        except:
                            minimum = 0         
            return minimum
        return 0
    
      
    def possible_reactions_backward(self,reaction_stoichiometry):   
        if str(reaction_stoichiometry) != "None":
            minimum = -1
            for metabolite in reaction_stoichiometry:
                
                #the quantity of backward reaction progressions is determined
                if metabolite[self.index_stoich_coef] > 0:   
                    if metabolite[self.index_compartment] == "c":
                        try:
                            temp = (self.cytoplasm_metabolites.loc[self.cytoplasm_metabolites["name"] == metabolite[self.index_name], "amount"]/abs(metabolite[self.index_stoich_coef])).iloc[0]
                            if temp < minimum:
                                minimum = temp
                        except:
                            minimum = 0          
                            
                    elif metabolite[self.index_compartment] == "m":
                        try:
                            temp = (self.cytoplasm_metabolites.loc[self.cytoplasm_metabolites["name"] == metabolite[self.index_name], "amount"]/abs(metabolite[self.index_stoich_coef])).iloc[0]
                            if temp < minimum:
                                minimum = temp
                        except:
                            minimum = 0         
            return minimum
        return 0
        
    
    def react(self,):   
        for index, row in self.reactions.iterrows():      
            enzyme_amount = enzyme_availability(row["Enzyme"])
            if str(row["Stoichiometry"]) != "None":
                if enzyme_amount > 0:              
                    extracellular_reaction = False  
    
                    #extracellular reactions are ignored
                    for metabolite in row["Stoichiometry"]:
                        if metabolite[self.index_compartment] != "c" and metabolite[self.index_compartment] != "m":  
                            extracellular_reaction = True
                 
                    if not extracellular_reaction:  
                        
                        #possible forward reactions with available enzymes are defined
                        if row["Direction"] == "Forward":
                            reaction_todo_forward = self.possible_reactions_forward(row["Stoichiometry"])*self.zeta_reaction_completeness   #!!! ticket 78
                            
                            if reaction_todo_forward > 0:    
                                if str(row["Forward kinetics"]) != "None":   
                                    enzyme_effectivness = float(self.enzymes.loc[row["Enzyme"], "effectiveness"])  
                                    max_enzyme_usage = float(reaction_todo_forward)*float(row["Forward kinetics"])*enzyme_effectivness  
                                    
                                    if reaction_todo_forward > max_enzyme_usage: 
                                        reaction_todo_forward = max_enzyme_usage
                                        
                                        if reaction_todo_forward > enzyme_amount:
                                            reaction_todo_forward = enzyme_amount
                                 
                                #calculates the concentration change after the reactions proceed
                                for metabolite in row["Stoichiometry"]:
                                    if metabolite[self.index_compartment] == "c":
                                        self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] += metabolite[self.index_stoich_coef]*reaction_todo_forward
                                    elif metabolite[self.index_compartment] == "m":
                                        self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] += metabolite[self.index_stoich_coef]*reaction_todo_forward
                                                
                            else:
                                print('ERROR: Incompatible kinetics and stoichiometry')
                        
                        #reaction progression is determined with reversible for which enzymes are available                  
                        elif row["Direction"] == "Reversible":
    
                            #reactions with potentially absorbable chemicals are evaluated
                            for index, row in self.reactions.iterrows():
                                if row["Stoichiometry"] != None:
                                    for chemical in row["Stoichiometry"]:
                                        if chemical[self.index_name] == potentially_absorbable_chemical["name"]:
                                            Q_numerator = 1
                                            Q_denominator = 1
                                            keq = 0
                                            current_Q = 0
                                            current_Q_keq = 0
                                            current_keq_displacement = 0
                                            
                                            if chemical[self.index_compartment] == "c":  
                                                if chemical[self.index_stoich_coef] > 0:  
                                                    Q_numerator *= pow(cytoplasm_metabolites.at[chemical[self.index_name], "amount"], chemical[self.index_stoich_coef])                          
                                                elif chemical[self.index_stoich_coef] < 0:    
                                                    Q_denominator *= pow(cytoplasm_metabolites.at[chemical[self.index_name], "amount"], -1*chemical[self.index_stoich_coef])                            
                    
                                            elif chemical[self.index_compartment] == "m":
                                                if chemical[self.index_stoich_coef] > 0:
                                                    Q_numerator *= pow(self.cytoplasm_metabolites.at[chemical[self.index_name], "amount"], chemical[self.index_stoich_coef])
                                                elif chemical[self.index_stoich_coef] < 0:
                                                    Q_denominator *= pow(self.cytoplasm_metabolites.at[chemical[self.index_name], "amount"], -1*chemical[self.index_stoich_coef])
                                                
                                            if Q_denominator != 0:
                                                current_Q = Q_numerator/Q_denominator
     
                                    
                                    
                            if Keq > current_Q:
                                
                                #forward reaction progression is considered
                                reaction_todo_forward = self.possible_reactions_forward(row["Stoichiometry"])*self.zeta_reaction_completeness 
                                
                                if reaction_todo_forward > 0:    
                                    if str(row["Forward kinetics"]) != "None":
                                        enzyme_effectivness = float(self.enzymes.loc[row["Enzyme"], "effectiveness"])
                                        max_enzyme_usage = float(reaction_todo_forward)*float(row["Forward kinetics"])*enzyme_effectivness 
                                        
                                        if reaction_todo_forward > max_enzyme_usage:
                                            reaction_todo_forward = max_enzyme_usage
                                            
                                            if reaction_todo_forward > enzyme_amount:
                                                reaction_todo_forward = enzyme_amount
                                                
                                    #calculates the concentration change after the reactions proceed     
                                    for metabolite in row["Stoichiometry"]:
                                        if metabolite[self.index_compartment] == "c":
                                            self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] += metabolite[self.index_stoich_coef]*reaction_todo_forward
                                        elif metabolite[self.index_compartment] == "m":
                                            self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] += metabolite[self.index_stoich_coef]*reaction_todo_forward
    
                            if Keq < current_Q:                       
    
                                #backward reaction progression is considered
                                reaction_todo_backward = self.possible_reactions_backward(row["Stoichiometry"])*self.zeta_reaction_completeness  
                                if reaction_todo_backward > 0:
                                    if str(row["Backward kinetics"]) != "None":
                                        enzyme_effectivness = float(self.enzymes.loc[row["Enzyme"], "effectiveness"])
                                        max_enzyme_usage = float(reaction_todo_backward)*float(row["Backward kinetics"])*enzyme_effectivness
                                        
                                        if reaction_todo_backward > max_enzyme_usage:
                                            reaction_todo_backward = max_enzyme_usage
                                            
                                            if reaction_todo_backward > enzyme_amount:
                                                reaction_todo_backward = enzyme_amount
                                                
                                    #calculates the concentration change after the reactions proceed                
                                    for metabolite in row["Stoichiometry"]:
                                        if metabolite[self.index_compartment] == "c":
                                            cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] -= metabolite[self.index_stoich_coef]*reaction_todo_backward
                                        elif metabolite[self.index_compartment] == "m":
                                            self.cytoplasm_metabolites.loc[metabolite[self.index_name], "amount"] -= metabolite[self.index_stoich_coef]*reaction_todo_backward
            
            