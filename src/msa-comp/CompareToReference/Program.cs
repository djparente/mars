/*
 * Created by SharpDevelop.
 * User: dparente
 * Date: 4/9/2012
 * Time: 15:48 PM
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Collections;
using System.Collections.Generic;
using DJPCommonBioinfo;

namespace CompareToReference
{
	class Program
	{
		public static void Main(string[] args)
		{
			string goldMSAPath = args[0];
			string candMSAPath = args[1];
			string prefix = "";//args[2];
			
			MSA goldMSA = new MSA(goldMSAPath);
			MSA candMSA = new MSA(candMSAPath);
			
			double score = compareMSAs(goldMSA, candMSA, prefix);
			
			Console.WriteLine(score);
		}

		
		//Compares two MSAs: one a gold standard, and the other ("candidate") an MSA 
		// produced as an attempt to re-create the gold standard.  
		public static double compareMSAs(MSA goldMSA, MSA candMSA, string realignedPrefix)
		{
			int hits = 0;
			int miss = 0;
			
			for(int i = 0; i < goldMSA.Rows; i++)
			{
				for(int j = 0; j < goldMSA.Rows; j++)	//Don't optimize this by setting j = i+1; an assumption relies on that not being that way
				{
					if( i == j )
						continue;
					
					//Console.Error.WriteLine("Comparing {0} and {1}", goldMSA.Names[i], goldMSA.Names[j]);
					
					//if( goldMSA.Names[i].StartsWith(realignedPrefix) && !goldMSA.Names[j].StartsWith(realignedPrefix) )
					if( i < j )
					{					
						
						//Console.Error.WriteLine("Comparing {0} and {1}", goldMSA.Names[i], goldMSA.Names[j]);
						int candI = getRowWithName(candMSA, goldMSA.Names[i]);
						int candJ = getRowWithName(candMSA, goldMSA.Names[j]);
						//Console.Error.WriteLine("Same as   {0} and {1}", candMSA.Names[candI], candMSA.Names[candJ]);
						
						
						
						List<Tuple<int, int>> goldMatch = getMatches(goldMSA, i, j);
						List<Tuple<int, int>> candMatch = getMatches(candMSA, candI, candJ);
						
						int curHit = 0;
						int curMiss = 0;
						
						matchScore(goldMatch, candMatch, out curHit, out curMiss);
						
						//Console.WriteLine("Result was : {0}", curMiss);	
						
						
						hits += curHit;
						miss += curMiss;
					}
				}
			}
			
			return (double)hits / (double)(hits+miss);
		}
		
		//Get matches betweent the specified rows
		public static List<Tuple<int, int>> getMatches(MSA msa, int r, int s)
		{
			int resR = 0;
			int resS = 0;
			
			List<Tuple<int,int>> rowMatches = new List<Tuple<int, int>>();
			
			for(int c = 0; c < msa.Cols; c++)
			{
				if( msa[r, c] != MSA.GAP )
					resR++;
				
				if( msa[s, c] != MSA.GAP )
					resS++;
				
				if( msa[s, c] != MSA.GAP && msa[r, c] != MSA.GAP )
					rowMatches.Add(new Tuple<int, int>(resR, resS));
					
			}
			
			return rowMatches;
		}
		
		public static void matchScore(List<Tuple<int, int>> goldMatches, List<Tuple<int, int>> candMatches, out int hit, out int miss)
		{
			hit = 0;
			miss = 0;
			
			foreach(Tuple<int, int> g in goldMatches )
			{
				if( candMatches.Contains(g) )
					hit++;
				else
					miss++;
			}
		}
		
		
		public static int getRowWithName(MSA msa, string name)
		{
			for(int i = 0; i < msa.Names.Length; i++)
				if( msa.Names[i] == name ) return i;
			
			throw new ArgumentException("Did not find row with that name: {0}", name);
		}
	}
}