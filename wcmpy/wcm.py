# -*- coding: utf-8 -*-
"""
@author: Andrew Freiburger & Ethan Chan
"""

#Import Statements
from scipy.constants import physical_constants, Boltzmann, R, femto, pico, nano
from math import pi, pow, sqrt, cos, radians, exp, floor
import codons
import pandas
import json, sys, os

def is_number(self,c):
    numbers = ".0123456789"
    if c in numbers:
        return True
    return False

def _sherical_radius(volume):
        return pow((volume*3)/(4*pi), 1/3)
    
def _spherical_surface_area(radius):
    return 4*pi*pow(radius, 2)

def _shell_volume(outer_radi, inner_radi):
        return (4*pi/3)*(outer_radi**3 - inner_radi**3) 

class BacMod():  # the code only currently supports spherical (-cocci) bacteria
    def __init__(self,):          
        self.codons = codons.Codons()
        self.bacterium = json.load(open(os.path.join(self.parameters['root_path'], 'parameters', 'bacteria', 'S_aureus.json')))
        self._primordial_conditions() 
        
        # import biochemial information
        self.amino_acid_synonyms = json.load(open(os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'amino_acid_synonyms.json')))
        self.enzyme_names = json.load(open(os.path.dirname(__file__), 'rosetta_stone', 'protein_names.json'))
        
        #environmental properties
        self.F = physical_constants['Faraday constant'] 
        self.time_log = 0  
        self.iteration = 0
        self.temperature = 298  #kelvins
        self.timestep_value = 10  #seconds   
        self.simulation_time = 1E2  #seconds
        self.simulation_timesteps = self.simulation_time / self.timestep_value  
        self.optimum_incubation_temperature = 310
        self.membrane_electrode_potential = 2


        self.cytoplasm_amino_acids = {"phenylalanine":0, "leucine":0, "isoleucine":0, "methionine":0, "valine":0, "serine":0, "proline":0, "threonine":0, "alanine":0, "tyrosine":0, "glutamine":0, "asparagine":0, "lysine":0, "aspartic acid":0, "glutamic acid":0, "glycine":0, "arginine":0, "histidine":0, "tryptophan":0, "cysteine":0}
        
        #calculation constants
        self.index_stoich_coef = 0
        self.index_compartment = 1
        self.index_name = 2
        self.zeta_reaction_completeness = 0.9   
        self.mitosis_time = 3869 #seconds
        self.s_phase_time = 15571 #seconds
        self.g_phase_time = 12960 #seconds
        
        # define initial DataFrames
        reactions = pandas.read_csv("./formatted_data/reactions.csv")
        reactions["Stoichiometry"] = reactions["Stoichiometry"].map(lambda stoichiometry: [[int(float(i.split(":")[self.index_stoich_coef])), i.split(":")[self.index_compartment], i.split(":")[self.index_name]] for i in str(stoichiometry).split(";")] if str(stoichiometry) != "None" and str(stoichiometry) != "nan" and str(stoichiometry) != "NaN" else None)
        
        self.cytoplasm_metabolites = pandas.read_csv("./formatted_data/molecules.csv" ) 
        #calculating the quantity of energetic molecules in the cytoplasm
        self.energetic_organism_molecules = []   
        for energetic_chemical in self.cytoplasm_metabolites["name"]:
            if self.cytoplasm_metabolites.loc[energetic_chemical, "energetic"] == "y":
                self.self.energetic_organism_molecules.append(energetic_chemical) 
        
        self.membrane_metabolites = self.bacterium['membrane_chemicals']
        self.genes = pandas.read_csv("./formatted_data/genes.csv")
        self.bulk_chemicals = pandas.DataFrame(columns = ["name", "concentration", "chemical_type", "barrier"])
        self.potentially_absorbed_chemicals = pandas.DataFrame(columns = ["name", "concentration", "chemical_type", "barrier"])
        self.enzymes = pandas.DataFrame(columns = ["name", "amount", "effectiveness"])
        self.antibiotics = pandas.DataFrame(columns = ["name", "concentration", "antibiotic_type", "magnitude", "target", "magnitude_lysis"])
        
    def _primordial_conditions(self,):
        # updating bacterial values
        membrane_thickness = 4*nano
        self.cell_volume = self.bacterium['cell_volume (fL)']*femto
        self.cell_radius = _sherical_radius(self.cell_volume)
        self.cell_surface_area = _spherical_surface_area(self.cell_radius)
        self.cell_absorption_volume = _shell_volume(self.cell_radius, (self.cell_radius-membrane_thickness))
        self.cell_mass = self.bacterium['cell_mass (pg)']*pico
        self.non_metabolic_mass_multiplier = self.bacterium['cell_mass (pg)']*pico / self.sum_metabolic_mass()
        
        #the extracellular concentrations are appended to the dataframe
        inital_bulk_chemicals = environment_harvested_data["bulk-chemicals"].split(",")
        for bulk_chemical in inital_bulk_chemicals:
            name, concentration, chemical_type, barrier = bulk_chemical.split(":")
            new_row = {"name":name, "concentration":float(concentration), "chemical_type":chemical_type, "barrier":float(barrier)}
            self.bulk_chemicals = bulk_chemicals.append(new_row, ignore_index = True)
    
        #the antibiotic concentrations are appended in the dataframe        
        inital_antibiotics = environment_harvested_data["antibiotics"].split(",")
        for antibiotic in inital_antibiotics:
            name, amount, antibiotic_type = antibiotic.split(":")
            new_row = {"name":name, "amount":float(amount), "antibiotic_type":antibiotic_type}
            self.antibiotics = antibiotics.append(new_row, ignore_index = True)
        
        #calculating the initial proportion of energetic chemicals
        primordial_quantity_molecules = 0
        primordial_energetic_quantity = 0
        for chemical in self.cytoplasm_metabolites["name"]:
            primordial_quantity_molecules += self.cytoplasm_metabolites.loc[chemical, "amount"]
            
            if chemical in self.energetic_organism_molecules:
                primordial_energetic_quantity += self.cytoplasm_metabolites.loc[chemical, "amount"] 
                
        self.primordial_energetic_proportion = primordial_energetic_quantity / primordial_quantity_molecules
        
        # default bacterial values
        self.cell_state = "primordial"
        self.cellular_vitality = "alive"
        self.timestep = 0
    
    def sum_metabolic_mass(self,):
        return (self.cytoplasm_metabolites["amount"] * self.cytoplasm_metabolites["Molecular weight"]).sum() + (self.cytoplasm_metabolites["amount"] * self.cytoplasm_metabolites["Molecular weight"]).sum()
    
    
    def average_metabolic_MW(self,):
        molecular_weights = 0
        molecules_total = 0
        for molecule in self.cytoplasm_metabolites:
            molecular_weights += self.cytoplasm_metabolites[molecule, 'Molecular weight']
            molecules_total += 1            
        return molecular_weights / molecules_total
    
    
    def cell_geometric_absorption(self,):
        #define the dimensions of the bacterium
        average_absorption_angle = 45    # degrees
        
        metabolic_mass = self.sum_metabolic_mass()
        self.cell_mass = metabolic_mass*self.non_metabolic_mass_multiplier
        self.cell_volume = self.cell_mass*((self.bacterium['cell_volume (fL)']*femto)/(self.bacterium['cell_mass (pg)']*pico))
        
        average_substrate_MW = self.average_metabolic_MW()
        root_mean_square_velocity = sqrt(3*Boltzmann*self.temperature/average_substrate_MW)
        absorption_shell_added_radius = root_mean_square_velocity * self.timestep_value  
                
        #calculate the cell surface area
        self.cell_radius = _sherical_radius(self.cell_volume)
        self.cell_surface_area = _spherical_surface_area(self.cell_radius)
    
        #calculate the shell volume of the potentially absorbable chemicals
        farthest_chemical_distance = self.cell_radius + absorption_shell_added_radius
        self.volume_absorption_shell = _shell_volume(farthest_chemical_distance, self.cell_radius)
        
        #calculate the proportion of absorbed chemicals within the shell volume
        component_height_velocty = absorption_shell_added_radius * cos(radians(average_absorption_angle))
        h = absorption_shell_added_radius-component_height_velocty
        absorption_cap_area = 2*absorption_shell_added_radius*h*pi
        possible_location_spherical_area = (4*pi*absorption_shell_added_radius**2)
        self.absorption_proportion = (absorption_cap_area / possible_location_spherical_area)
        
    def _reset(self,):
        self.enzymes["effectiveness"] = 1.0
        self.recycle_amino_acids()
        self.cell_geometric_absorption()
                 
    def _setup(self,):
        enzymes["effectiveness"] = 1.0
        self.cell_geometric_absorption()
        
        #the non-metabolic mass at the initial time calculated as a constant
        if self.iteration == 0:
            metabolic_mass = self.sum_metabolic_mass()
            self.non_metabolic_mass_multiplier = self.bacterium['cell_mass (pg)']*pico / metabolic_mass  
    
    
    def antibiotic_effects(self,):
        for row in self.antibiotics.itertuples(self,):
            row_index = row.Index     
            antibiotic_type = row.antibiotic_type
            antibiotic_amount = row.amount
            antibiotic_name = row.name
            antibiotic_magnitude = row.magnitude
            antibiotic_target = row.target
            
            if antibiotic_type == "inhibit":
                self.enzymes.at[antibiotic_target, "effectiveness"] -= antibiotic_amount*antibiotic_magnitude
                if self.enzymes.at[antibiotic_target, "effectiveness"] < 0:
                    antibiotic_target.at[antibiotic_target, "effectiveness"] = 0
        
    

    
    
    def binary_fission(self,):        
        cell_cycle_duration = self.mitosis_time + self.s_phase_time + self.g_phase_time
     
        replicated_times = floor(self.time_log/cell_cycle_duration)
        
        if self.iteration == 0:
            print("G_phase: " + (self.time_log - replicated_times*cell_cycle_duration) + " seconds")
        
        if self.time_log - cell_cycle_duration * replicated_times == self.g_phase_time:  
            print("S_phase: " + (self.time_log - replicated_times*cell_cycle_duration) + " seconds")  
        
        elif self.time_log - cell_cycle_duration * replicated_times > self.g_phase_time and self.time_log == self.s_phase_time + self.g_phase_time: 
            print("Mitosis: " + (self.time_log - replicated_times*cell_cycle_duration) + " seconds")
        
        elif self.time_log - cell_cycle_duration * replicated_times > self.g_phase_time + self.s_phase_time and self.time_log == cell_cycle_duration:  
            print("replication: %d" % (replicated_times))
            self.cell_state = "primordial"   
            
     
    def _death(self,):  
        #calculate the region of membrane oxidation
        oxidation_angle = 5
        component_height = self.cell_radius * cos(radians(oxidation_angle))
        h = self.cell_radius-component_height
        oxidized_cap_area = 2*self.cell_radius*h*pi
        oxidized_area_proportion = (oxidized_cap_area / self.cell_surface_area)
        
        #determine the quantities of chemicals in the membrane compartment
        oxidized_membrane_chemicals = 0
        unsaturated_membrane_chemicals = 0    
        
        for index, metabolite in self.cytoplasm_metabolites.iterrows():
            if metabolite["unsaturated"]:   
                unsaturated_membrane_chemicals += metabolite["amount"]
                
            if metabolite["oxidized"]:
                oxidized_membrane_chemicals += metabolite["amount"]
                   
        #determine fatality from lysis
        oxidation_vitality = (oxidized_membrane_chemicals / unsaturated_membrane_chemicals)/oxidized_area_proportion
        if oxidation_vitality >= 1:
            cellular_vitality = "dead"
            
        #determine fatality from wasting 
        wasting_death_proportion = 1/3
        if self.cell_mass <= self.bacterium['cell_mass (pg)']*pico*wasting_death_proportion:
            cellular_vitality = "dead"
            
        #determine fatality from energetic death
        quantity_energetic_molecules = 0
        total_quantity_molecules = 0
        
        for chemical in self.cytoplasm_metabolites["name"]:
            total_quantity_molecules += self.cytoplasm_metabolites[chemical]["amount"]
            
            if chemical in self.energetic_organism_molecules:
                quantity_energetic_molecules += self.cytoplasm_metabolites[chemical]["amount"]        
                
        proportion_energetic_molecules = quantity_energetic_molecules / total_quantity_molecules
        
        if proportion_energetic_molecules < float(self.primordial_energetic_proportion/2) :
            cellular_vitality = "dead"
            
        #simulation end from death
        if cellular_vitality == "dead": 
            print ("The bacterium has died.\nThe simulation is over.")  
            sys.exit()    
    
    def translation(self,):    
        for index, gene in self.genes.iterrows():
            quantity_of_codons = len(gene) / 3
            time_per_codon = self.mitosis_time / quantity_of_codons   
            codons_per_timestep = self.timestep_value / time_per_codon 
            
            # the created amino acids in each timestep are determined
            amino_acid_sequence = self.codons(gene)
            created_amino_acids = amino_acid_sequence[self.timestep*codons_per_timestep:(self.timestep+1)*codons_per_timestep]
            for amino_acid in created_amino_acids:
                amino_acid_name = self.amino_acid_synonyms[amino_acid]['name']
                
                if self.cytoplasm_amino_acids[amino_acid_name] >= 1:  
                    self.cytoplasm_amino_acids[amino_acid_name] -= 1
                
            #the end of sequences in determined and the synthesized protein is appended to the list of enzymes
            if len(gene) <= codons_per_timestep*self.timestep:
                enzyme_name = self.enzyme_names[amino_acid_sequence]
                self.enzymes.append({"name":enzyme_name, "amount":1, "effectiveness":1.0} , ignore_index = True) 

                #the genes with numerical names are reformatted                     
                stripped_gene = "".join(filter(is_number, str(enzyme_name))).lstrip("0")
                
                for index, enzyme in self.enzymes.iterrows():
                    if "".join(filter(is_number, str(enzyme["name"]))).lstrip("0") == stripped_gene:    
                        self.enzymes.loc[index, "amount"] += 1 
    
    def recycle_amino_acids(self,):
        #every codon is read and the corresponding amino acids are added to the cytoplasmic concentrations
        enzymes_halflife = 5000    #seconds
        degradation_time = self.time_log / 2
        for enzyme_index, enzyme in self.enzymes.iterrows():
    
            #the decayed enzymes are subtracted from total quantity of enzymes
            enzymes_to_recycle = enzyme["amount"]*pow(1/2, degradation_time/enzymes_halflife) 
            enzyme["amount"] -= enzymes_to_recycle
            
            if enzymes_to_recycle > 0:  
                gene_sequence = self.genes[self.genes["name"] == enzyme["name"], "Sequence"] 
    
                for i in range(int(len(gene_sequence)/3)): 
                     codon = gene_sequence[i*3:(i*3)+3]
                     amino_acid = self.classify_codon(codon)
                     
                     if amino_acid == "break":
                         break
                     else:
                         self.cytoplasm_amino_acids[amino_acid] += 1
        

    def enzyme_availability(self,enzyme_name):
        if not enzyme_name == "None":
            try:
                return self.enzymes.loc[self.enzymes["name"] == enzyme_name, "amount"].iloc[0]
            except:
                return 0
    
        else:
            return -1
            
      
    
    def state(self,):
        self.death()
        self.binary_fission()
        
    def plotting(self,):
        pass
    
    def export(self,):
        pass
    
    
            
    def simulate(self,):    
        self.format_dataframes()
        self.import_environment(environment_path)
        self.setup()
        
        #Simulation starting point, loop through given iterations
        for self.timestep in range(int(self.simulation_timesteps)):
            self.initial_conditions()
            self.absorption()
            self.translation()
            self.antibiotic_effects()
            self.react()
            self.state()
            self._reset()
            self.iteration += 1
            self.time_log += self.timestep_value
            
            if self.iteration == 1:
                self.cell_state = "adult"
            
            if self.iteration == self.simulation_timesteps:  #!!!
                print("The simulation is over.")
                self.plotting()  #!!!
            
            
            
# execute the simulation
bacterial_model = BacMod()
bacterial_model.main()