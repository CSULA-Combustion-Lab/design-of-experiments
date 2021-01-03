# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 12:10:13 2019

@author: boss
"""

import cantera
import matplotlib.pyplot as plt 
import numpy as np
import math
import os
import datetime
import time
#import pickle
#import cProfile
#from statistics import mode
cantera.suppress_thermo_warnings()


def plot_directory():
    """ Creates directories inside the working directory for all the plots in 
        the script. current_directory/Plots/date/file_number/plot_type/file """ 
        
    file_num = 1
    current_dir = os.getcwd()
    ranks_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Ranking')
    while os.path.exists(ranks_dir) == True: 
        file_num = file_num+1
        ranks_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Ranking')    
    if not os.path.exists(ranks_dir):
        os.makedirs(ranks_dir)
    score2_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Score2')
    if not os.path.exists(score2_dir):
        os.makedirs(score2_dir)   
    score_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Score')
    if not os.path.exists(score_dir):
        os.makedirs(score_dir) 
    data_dir = os.path.join(current_dir,'Plots/'+date+'/'+str(file_num)+'/Data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)        
    return file_num


def interested_rxn_ranking(all_ranks):
    """ Takes the rankings of all reactions and returns a list of the ranking 
    of the reactions of interest for each species for each condition and the 
    initial inputs for that condition""" 
    all_ranks = np.array(all_ranks)
    indices = rxns_interested_indices(rxns,interested_rxns) 
    time_steps = len(senstime)
    species_sens = [None]*time_steps
    rank_across_spec =[None]*len(indices)
    ranks_per_cond = []
    rank = [None]*5
    rank[0] = phi
    rank[1] = methane
    rank[2] = temp
    rank[3] = pressure*101.325
    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(indices)):
            sense = all_ranks[:,indices[j]] 
            for k in range(0,time_steps):
                species_sens[k] = sense[i+k*len(SpecificSpecies)]
            rank_across_spec[j] = min(species_sens)
        ranks_per_cond.append([x for x in rank_across_spec])
        rank_plot.append([x for x in rank_across_spec])
    rank[4] = np.array(ranks_per_cond)
    rank_score.append([x for x in rank])
    return rank_score, rank_plot


def rxns_interested_indices(rxns,interested_rxns):
    """ Creates an array of indicies of the reactions of interest with respect
    to the list of all reactions"""
    rxns = rxns.tolist()
    rxn_indx = []
    for rxn in interested_rxns:
        ind = rxns.index(rxn)
        rxn_indx.append(ind)
    return rxn_indx

def ranking_plots(rank_score,rank_plot):
    
    Phi_ary=[None]*len(rank_score)
    CH4_ary=[None]*len(rank_score)
    T_ary=[None]*len(rank_score)
    P_ary=[None]*len(rank_score)
    for cond in range(0,len(rank_score)):    
        Phi_ary[cond]=rank_score[cond][0]
        CH4_ary[cond]=rank_score[cond][1]
        T_ary[cond]=rank_score[cond][2]
        P_ary[cond]=rank_score[cond][3]    
    Phi_ary=np.array(Phi_ary)
    CH4_ary=np.array(CH4_ary)
    T_ary=np.array(T_ary)
    P_ary=np.array(P_ary)
    condition_num = int(len(rank_plot)/len(SpecificSpecies))
    for rxn in interested_rxns:
        ranks = [None]*condition_num
        fig1 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax1 = fig1.add_subplot(111)
        ax1.grid()
        ax1.set_ylabel('Ranking')
        ax1.set_xlabel('Pressure [kPa]')
        ax1.set_title('Rxn '+str(rxn)+' Ranking vs Presure')
        fig2 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax2 = fig2.add_subplot(111)
        ax2.grid()
        ax2.set_ylabel('Ranking')
        ax2.set_xlabel('Temperature [K]')
        ax2.set_title('Rxn '+str(rxn)+' Ranking vs Temperature')
        fig3 = plt.figure()
        plt.ylim(len(rxns)+1,0)
        ax3 = fig3.add_subplot(111)
        ax3.grid()
        ax3.set_ylabel('Ranking')
        ax3.set_xlabel('Methane.X')
        ax3.set_title('Rxn '+str(rxn)+' Ranking vs Methane Mole Fraction')
        fig4 = plt.figure()
        plt.ylim(len(rxns)+1,0) 
        ax4 = fig4.add_subplot(111)
        ax4.grid()
        ax4.set_ylabel('Ranking')
        ax4.set_xlabel('Equivalence Ratio')
        ax4.set_title('Rxn '+str(rxn)+' Ranking vs Equivalence Ratio')
        marker=['o','s','v','^','>']
        color=['b','darkorange','g','r','blueviolet']
        for species, mrk, c in zip(SpecificSpecies,marker,color):
            for i in range(0,condition_num):
                ranks[i] = rank_plot[i*len(SpecificSpecies)+
                                       SpecificSpecies.index(species)][interested_rxns.index(rxn)]
            ranks = np.array(ranks)
            ax1.scatter(P_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig1.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig1.tight_layout()
            fig1.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Pressure.png',bbox_inches = 'tight')
            ax2.scatter(T_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig2.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig2.tight_layout()
            fig2.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Tempurature.png',bbox_inches = 'tight')
            ax3.scatter(CH4_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig3.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig3.tight_layout()
            fig3.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Methane%.png',bbox_inches = 'tight')
            ax4.scatter(Phi_ary,ranks,label=species,marker=mrk,facecolor='none'
                        ,edgecolor=c)
            fig4.legend(bbox_to_anchor=(0.95, 1), loc=2)
            fig4.tight_layout()
            fig4.savefig('Plots/'+date+'/'+str(file_num)+'/Ranking/Rxn'+str(rxn)+'_Phi.png',bbox_inches = 'tight')
            
if __name__ == "__main__":
    ranking_plots(rank_score,rank_plot)
