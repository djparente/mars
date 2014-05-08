/*
 * Daniel J Parente
 * Swint-Kruse Laboratory
 * University of Kansas Medical Center
 *
 * Date: 12/15/2011
 * Time: 14:26 PM
 * 
 */
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DJPCommonBioinfo;

namespace ConstrainedNeedle
{
	enum TraceSource { Gap, Insert, Match, Done };
	struct TraceCell
	{
		private double score;
		private TraceSource source;
		
		//Protect private fields and make object immutable
		public double Score { get { return this.score; } }
		public TraceSource Source { get { return this.source; } }
		
		public TraceCell(double score, TraceSource source)
		{
			this.score = score;
			this.source = source;
		}
	}
	
	class Program
	{
		public static void Main(string[] args)
		{
			try
			{
				runProgram(args);
			}
			catch(Exception e)
			{
				displayError("An unexpected error occured: " + e.Message);
				Environment.Exit(-1);
			}
		}
		
		private static void runProgram(string[] args)
		{
			if( args.Length != 3 )
			{
				displayError("Expected three arguments, but was given " + args.Length + " arguments instead");
				usage();
				Environment.Exit(-254);
			}
			
			Console.Error.WriteLine("Parsing reference MSA");
			MSA refer = new MSA(args[0]);
			
			Console.Error.WriteLine("Parsing query MSA");
			MSA query = new MSA(args[1]);
			
			Console.Error.WriteLine("Parsing similarity matrix");
			SimilarityMatrix simMat = new SimilarityMatrix(args[2]);
			
			int minNameLength = Math.Max(refer.getLongestNameLength(), query.getLongestNameLength());
			
			List<Tuple<int, int>> seqPairs = findEqualSequences(refer, query);		
			//List<Tuple<int, int>> seqPairs = findEqualSequencesWithTree(refer, query);
			
			Console.Error.WriteLine("Constraints:");
			for(int i = 0; i < seqPairs.Count; i++)
				Console.Error.WriteLine("+Constraint {0} <---> {1}", refer.Names[seqPairs[i].Item1],query.Names[seqPairs[i].Item2]);
			
			Console.Error.WriteLine("Found a total of {0} constraints", seqPairs.Count);
			
			if( seqPairs.Count == 0 )
			{
				displayError("WARNING: NO CONSTRAINTS DETECTED!");
				//Console.Error.WriteLine("WARNING: NO CONSTRAINTS DETECTED!");
			}
			
//			refer.displayConsole(minNameLength);
//			query.displayConsole(minNameLength);
			
			Console.Error.WriteLine("Building constraint matrix");
			double[,] constraint = GetConstraintMatrix(refer, query, seqPairs);
			Console.Error.WriteLine("Building similarity matrix");
			//double[,] similarity = getSimilarityConstraintMatrix(refer, query, simMat);
			double[,] similarity = getSimilarityConstraintMatrixByFreq(refer, query, simMat, true);
			
			Console.Error.WriteLine("Merging into single energy matrix");
			scaleOffsetConstraintMatrix(similarity, .1/(simMat.Max - simMat.Min), -simMat.Min);
			constraint = sumConstraints(constraint, similarity);
			
//			for(int i = 0; i < constraint.GetLength(0); i++)
//			{
//				for(int j = 0; j < constraint.GetLength(1); j++)
//				{
//					Console.Write(constraint[i,j].ToString("F2") + "\t");
//				}
//				Console.WriteLine();
//			}
			
			Console.Error.WriteLine("Performing NW-alignment");
			List<TraceSource> result = needleAlignSequential(constraint, 0.01, 0.01);//Used to be 0.01, 0.01 during normal operation 
//			foreach(TraceSource source in result )
//				Console.WriteLine(source.ToString());
			
			//MSA newRefer = applyAlignmentToMSA(result, refer, false);
			//MSA newQuery = applyAlignmentToMSA(result, query, true);
			
			
			
			//newRefer.displayConsole(minNameLength);
			//newQuery.displayConsole(minNameLength);
			//Console.WriteLine(newRefer.getFastaString());
			//Console.WriteLine(newQuery.getFastaString());
			Console.Error.WriteLine("Building final MSA");
			MSA final = buildFinalMSA(result, refer, query, seqPairs, true);
			Console.WriteLine(final.getFastaString());
		}
		
		public static void usage()
		{
			Console.Error.WriteLine("Usage: ./MARS.exe [Query alignment path] [Reference Alignment Path] [Similarity Matrix Path] > [Output file]");
		}
		
		public static void displayError(string message)
		{
			ConsoleColor prevColor = Console.ForegroundColor;
			Console.ForegroundColor = ConsoleColor.Red;
			
			Console.Error.WriteLine(message);
			
			Console.ForegroundColor = prevColor;
		}
		
		//Gap and insert penalty should be positive numbers (which will be SUBTRACTED from quality score)
		public static List<TraceSource> needleAlignSequential(double[,] scoreMatrix, double gapPenalty, double insertPenalty)
		{
			int rows = scoreMatrix.GetLength(0)+1;
			int cols = scoreMatrix.GetLength(1)+1;
			
			TraceCell[,] traceMatrix = new TraceCell[rows, cols];
			//Initialize the trace matrix
			traceMatrix[0,0] = new TraceCell(0, TraceSource.Done);
			for(int r = 1; r < rows; r++)
				traceMatrix[r, 0] = new TraceCell(traceMatrix[r-1 ,0].Score - insertPenalty, TraceSource.Insert);
			for(int c = 1; c < cols; c++)
				traceMatrix[0, c] = new TraceCell(traceMatrix[0, c-1].Score - gapPenalty, TraceSource.Gap);
			
//			//Calculate entries diagonally
//			for(int i = 0; i < rows + cols - 1; i++)
//			{				
//				int kMin = ( i <= cols ) ? 1 : 0;
//				
//				int initR = Math.Max(0, i - cols + 1);
//				int initC = Math.Min(i, cols - 1);
//				
//				int kMax = Math.Max(0, Math.Min(initC - 1, rows - initR));
//				
//				Parallel.For(kMin, kMax, k =>
//				//for(int k = kMin; k < kMax; k++)
//				{
//					int r = initR + k;
//					int c = initC - k;
					
			for(int r = 1; r < rows; r++)
			{
				for(int c = 1; c < cols; c++)
				{
					double matchScore = traceMatrix[r-1, c-1].Score + scoreMatrix[r-1,c-1];
					double gapScore = traceMatrix[r, c-1].Score - gapPenalty;
					double insertScore = traceMatrix[r-1, c].Score - insertPenalty;
					
					double maxScore = Math.Max(matchScore, Math.Max(gapScore, insertScore));
					
					//Prefer in this order: match, gap, insert
					if( maxScore == matchScore )
						traceMatrix[r,c] = new TraceCell(matchScore, TraceSource.Match);
					else if( maxScore == gapScore )
						traceMatrix[r,c] = new TraceCell(gapScore, TraceSource.Gap);
					else
						traceMatrix[r,c] = new TraceCell(insertScore, TraceSource.Insert);
				}
			}
			
//			for(int r = 0; r < rows; r++)
//			{
//				for(int c = 0; c < cols; c++)
//					Console.Write(traceMatrix[r,c].Score.ToString("F2").PadLeft(5,' ') + "\t");
//				Console.WriteLine();
//			}
			
			//Map the pathway
			List<TraceSource> result = new List<TraceSource>();
			int curR = rows-1;
			int curC = cols-1;
			while( curR > 0 || curC > 0 )
			{
				result.Add(traceMatrix[curR, curC].Source);
				
				switch( traceMatrix[curR, curC].Source )
				{
					case TraceSource.Match:
						curR -= 1;
						curC -=1;
						break;
					case TraceSource.Gap:
						curC -=1;
						break;
					case TraceSource.Insert:
						curR -= 1;
						break;
					case TraceSource.Done:
						throw new InvalidOperationException();
				}
					
			}
			
			//Reverse the array
			result.Reverse();
			return result;
		}
	
		#region Constraint-extraction methods
		public static double[,] GetConstraintMatrix(MSA refMSA, MSA queryMSA, List<Tuple<int, int>> seqPairs)
		{
			//Find equivalent sequences
			//List<Tuple<int, int>> seqPairs = findEqualSequences(refMSA, queryMSA);
			
			//Set up for constraint determination
			List<int>[] colConstraints = new List<int>[queryMSA.Cols];
			for(int i = 0; i < colConstraints.Length; i++)
				colConstraints[i] = new List<int>();
			
			//Actually calculate the constraints as a list
			foreach(Tuple<int, int> pair in seqPairs)
				addColConstraint(refMSA, queryMSA, pair.Item1, pair.Item2, colConstraints);
			
			//Generate a constraint matrix
			double[,] constraintMatrix = new double[queryMSA.Cols, refMSA.Cols];
			for(int i = 0; i < colConstraints.Length; i++)
			{
				for(int j = 0; j < colConstraints[i].Count; j++)
					constraintMatrix[i, colConstraints[i][j]]++;
			}
			
			return constraintMatrix;
		}
		
		public static double[,] getSimilarityConstraintMatrix(MSA refMSA, MSA queryMSA, SimilarityMatrix simMatrix)
		{
			double[,] constraintMatrix = new double[queryMSA.Cols, refMSA.Cols];
			for(int rCol = 0; rCol < refMSA.Cols; rCol++)
			{
				Console.Error.WriteLine("Calculating row {0} of {1}", rCol+1, refMSA.Cols );
				for(int qCol = 0; qCol < queryMSA.Cols; qCol++)
					constraintMatrix[qCol, rCol] = columnwiseSimilarity(rCol, qCol, refMSA, queryMSA, simMatrix, true, 0);
			}
			
			return constraintMatrix;
		}
		
		public static double[,] getSimilarityConstraintMatrixByFreq(MSA refMSA, MSA queryMSA, SimilarityMatrix simMatrix, bool ignoreGaps)
		{
			Dictionary<char, double>[] qFreqs = new Dictionary<char, double>[queryMSA.Cols];
			
			double[,] constraintMatrix = new double[queryMSA.Cols, refMSA.Cols];
			for(int rCol = 0; rCol < refMSA.Cols; rCol++)
			{
				//Console.Error.WriteLine("Calculating row {0} of {1}", rCol+1, refMSA.Cols );
				Dictionary<char, double> rFreq = getColFreq(refMSA, rCol, ignoreGaps);
				for(int qCol = 0; qCol < queryMSA.Cols; qCol++)
				{
					if( qFreqs[qCol] == null)
						qFreqs[qCol] = getColFreq(queryMSA, qCol, ignoreGaps);
					
					constraintMatrix[qCol, rCol] = compareColsByFreq(rFreq, qFreqs[qCol], simMatrix);
				}
			}
			
			return constraintMatrix;
		}
		
		public static Dictionary<char, double> getColFreq(MSA msa, int col, bool ignoreGaps)
		{
			int count = 0;
			Dictionary<char, double> result = new Dictionary<char, double>();
			
			List<char> chars = new List<char>();
			for(int i = 0; i < msa.Rows; i++)
			{
				if( ignoreGaps && msa[i,col] == MSA.GAP )
					continue;
				
				if( !result.ContainsKey(msa[i, col]) )
				{
					result.Add(msa[i,col], 1);
					chars.Add(msa[i,col]);
				}
				else
					result[msa[i,col]]++;
				
				count++;
			}
			
			foreach( char c in chars )
				result[c] /= count;
				
			
			return result;
		}
		
		public static double compareColsByFreq(Dictionary<char, double> col1, Dictionary<char, double> col2, SimilarityMatrix simMatrix)
		{
			double result = 0;
			foreach(char c1 in col1.Keys)
				foreach(char c2 in col2.Keys)
					result += simMatrix[c1,c2]*col1[c1]*col2[c2];
			
			return result;
		}
		
		public static double columnwiseSimilarity(int rCol, int qCol, MSA rMSA, MSA qMSA, SimilarityMatrix simMatrix, bool ignoreGaps, double NO_RESULT)
		{
			double result = 0;
			int count = 0;
			for(int rRow = 0; rRow < rMSA.Rows; rRow++)
			{
				if( ignoreGaps && rMSA[rRow,rCol] == MSA.GAP )
					continue;
				
				for(int qRow = 0; qRow < qMSA.Rows; qRow++)
				{			
					if( ignoreGaps && qMSA[qRow,qCol] == MSA.GAP )
						continue;
					
					result += simMatrix[rMSA[rRow,rCol],qMSA[qRow,qCol]];
					count++;
				}				
			}
			
			if( count == 0 )
				return NO_RESULT;
			
			result /= count;
			
			return result;
		}
		
		//Offset first, then scale
		public static void scaleOffsetConstraintMatrix(double[,] constraintMatrix, double scaleFactor, double offset)
		{
			int rows = constraintMatrix.GetLength(0);
			int cols = constraintMatrix.GetLength(1);
			
			for(int r = 0; r < rows; r++)
			{
				for(int c = 0; c < cols; c++)
				{
					constraintMatrix[r,c] = (constraintMatrix[r,c] + offset)*scaleFactor;
				}
			}
		}
		
		public static double[,] sumConstraints(double[,] c1, double[,] c2)
		{
			int c1Row = c1.GetLength(0);
			int c1Col = c1.GetLength(1);
			
			int c2Row = c1.GetLength(0);
			int c2Col = c2.GetLength(1);
			
			if( c1Row != c2Row || c1Col != c2Col )
				throw new ArgumentException("Constraint matrices must be of the same size");
			
			double[,] result = new double[c1Row, c1Col];
			for(int r = 0; r < c1Row; r++)
			{
				for(int c = 0; c < c2Col; c++)
				{
					result[r,c] = c1[r,c] + c2[r,c];
				}
			}
			
			return result;
		}
				
		#region Deprecated computePercentages method
		/*
		//Returns two arrays (with related data at matched indices)
		// target[i] specifies the column in MSA1 (ref) that corrosponds to position i in MSA2 (query)
		// percentTargetCoverage[i] percent of the time.
		static void computePercentages(List<int> constraints, out List<int> target, out List<int> count, out List<double> percentTargetCoverage)
		{
			constraints.Sort();
			
			//Set up the return variables
			target = new List<int>();
			count = new List<int>();
			percentTargetCoverage = new List<double>();
			
			//All processing done if no constraints
			if( constraints.Count == 0 )
				return;
			
			//Initialize the targetList and countList 
			count.Add(1);
			target.Add(constraints[0]);
			int numDistinctTargets = 1;
			for(int i = 1; i < constraints.Count; i++)
			{
				if( constraints[i] == constraints[i-1] )
					count[numDistinctTargets-1]++;
				else
				{
					numDistinctTargets++;
					target.Add(constraints[i]);
					count.Add(1);
				}
			}
			
			//Calculate percents based on the counts
			for(int i = 0; i < target.Count; i++)
				percentTargetCoverage.Add((double)(count[i]) / constraints.Count);
		}
		*/
		#endregion
		
		//Adds a constraint to the colHomology list for all residues along a matched
		// sequence
		static void addColConstraint(MSA m1, MSA m2, int r1, int r2, List<int>[] colHomology)
		{			
			int colPtr1 = -1;
			int colPtr2 = -1;
			
			nextResidue(m1, r1, ref colPtr1);
			nextResidue(m2, r2, ref colPtr2);
			
			while( colPtr1 < m1.Cols && colPtr2 < m2.Cols )
			{
				colHomology[colPtr2].Add(colPtr1);
				nextResidue(m1, r1, ref colPtr1);
				nextResidue(m2, r2, ref colPtr2);
			}
		}
		
		//Increment the pointer to the next residue in the specified row of the MSA
		static void nextResidue(MSA msa, int row, ref int colPtr)
		{
			colPtr++;
			while( colPtr < msa.Cols && msa[row, colPtr] == '-' )
				colPtr++;
		}
		
		//Returns a list of sequences which are equivalent
		static List<Tuple<int, int>> findEqualSequences(MSA m1, MSA m2)
		{
			List<Tuple<int, int>> seqPairs = new List<Tuple<int, int>>();
			
			Dictionary<int, string> r2Cache = new Dictionary<int, string>(m2.Rows);
			
			for(int r1 = 0; r1 < m1.Rows; r1++)
			{
				string r1Seq = m1.extractUngappedSequence(r1);
				for(int r2 = 0; r2 < m2.Rows; r2++)
				{
					//Get the MSA2 sequence from the cache (or compute it if first encounter)
					string r2Seq;
					if( !r2Cache.TryGetValue(r2, out r2Seq) )
					{
						r2Seq = m2.extractUngappedSequence(r2);
						r2Cache.Add(r2, r2Seq);
					}
					
					//Compare sequences
					if( r1Seq == r2Seq )
						seqPairs.Add(new Tuple<int, int>(r1, r2));
				}
				
			}
			
			return seqPairs;
		}
		
//		public static List<Tuple<int, int>> findEqualSequencesWithTree(MSA m1, MSA m2)
//		{
//			List<Tuple<int, int>> seqPairs = new List<Tuple<int, int>>();
//			ContainsTree cTree = new ContainsTree(m1);
//			
//			for(int i = 0; i < m2.Rows; i++)
//			{
//				List<int> cRow = cTree.contains(m2, i);
//				if( cRow == null )
//					continue;
//				
//				for(int j = 0; j < cRow.Count; j++)
//					seqPairs.Add(new Tuple<int, int>(j, i));
//			}
//			
//			return seqPairs;
//		}
		#endregion
		
		public static MSA applyAlignmentToMSA(List<TraceSource> alignment, MSA msa, bool gapIsGap)
		{
			TraceSource gapCmd = gapIsGap ? TraceSource.Gap : TraceSource.Insert;
			
			int newGaps = 0;
			foreach( TraceSource cmd in alignment )
				if( cmd == gapCmd )
					newGaps++;
			
			int newRows = msa.Rows;
			int newCols = msa.Cols + newGaps;
			
			char[,] newAln = new char[newRows, newCols];
			int usedPtr = 0;
			bool alwaysGap;
			for(int i = 0; i < alignment.Count; i++)
			{
				alwaysGap = alignment[i] == gapCmd;
				for(int r = 0; r < newRows; r++)
						if( alwaysGap )
							newAln[r,i] = MSA.GAP;
						else
							newAln[r,i] = msa[r,usedPtr];
							
				
				if( !alwaysGap )
					usedPtr++;
			}
			
			return new MSA(newAln, msa.Names);
		}
		
		public static MSA buildFinalMSA(List<TraceSource> alignment, MSA m1, MSA m2, List<Tuple<int, int>> equivSeqs, bool ignoreDuplicates)
		{
			bool[] ignoreSeq = new bool[m2.Rows];
			int numSeqIgnore = 0;
			if( ignoreDuplicates )
			{
				for(int i = 0; i < equivSeqs.Count; i++)
				{
					if( ignoreSeq[equivSeqs[i].Item2] != true )
					{
						numSeqIgnore++;
						ignoreSeq[equivSeqs[i].Item2] = true;
					}
				}
			}
			
			int newRows = m1.Rows + m2.Rows - numSeqIgnore;
			int newCols = alignment.Count;
			
			char[,] newAln = new char[newRows, newCols];
			
			//Add MSA1
			int usedPtr = 0;
			bool alwaysGap;
			for(int i = 0; i < alignment.Count; i++)
			{
				alwaysGap = alignment[i] == TraceSource.Insert;
				for(int r = 0; r < m1.Rows; r++)
						if( alwaysGap )
							newAln[r,i] = MSA.GAP;
						else
							newAln[r,i] = m1[r,usedPtr];
							
				
				if( !alwaysGap )
					usedPtr++;
			}
			
			//ADD MSA2
			usedPtr = 0;
			for(int i = 0; i < alignment.Count; i++)
			{
				alwaysGap = alignment[i] == TraceSource.Gap;
				int curRow = m1.Rows;
				for(int r = 0; r < m2.Rows; r++)
				{
					if( ignoreSeq[r] )
						continue;
					
					if( alwaysGap )
						newAln[curRow,i] = MSA.GAP;
					else
						newAln[curRow,i] = m2[r,usedPtr];
					
					curRow++;
				}
							
				if( !alwaysGap )
					usedPtr++;
			}
			
			//Construct the names list
			string[] newNames = new string[newRows];
			Array.Copy(m1.Names, newNames, m1.Names.Length);
			
			int rowPtr = m1.Rows;
			for(int r = 0; r < m2.Rows; r++)
			{
				if( ignoreSeq[r] )
					continue;
				
				newNames[rowPtr] = m2.Names[r];
				
				rowPtr++;
			}
			
			return new MSA(newAln, newNames);
		}
	}
	
	
}