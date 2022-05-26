package cn.edu.jxufe.agl;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Random;
 
public class GeneticAlgorithmTSP {
 
    /**
     * ��Ⱥ��ģ
     */
	private int scale;
	
	/**
	 * ����������Ⱦɫ�峤��
	 */
	private int cityNum;
	
	/**
	 * ��ֳ����
	 */
	private int MAX_GEN;
	
	/**
	 * ���о������
	 */
	private int[][] distances;
	
	/**
	 * ��ѳ��ִ���
	 */
	private int bestT;
	
	/**
	 * ��ѳ���
	 */
	private int bestLength;
	
	/**
	 * ��Ѹ��壨·����
	 */
	private int[] bestTour;
 
	/**
	 * ��ʼ��Ⱥ��������Ⱥ��������ʾ��Ⱥ��ģ��һ�д���һ�����壬��Ⱦɫ�壬�б�ʾȾɫ�����Ƭ��
	 */
	private int[][] oldPopulation;
	
	/**
	 * �µ���Ⱥ���Ӵ���Ⱥ
	 */
	private int[][] newPopulation;
	
	/**
	 * ��Ⱥ��Ӧ�ȣ���ʾ��Ⱥ�и����������Ӧ��
	 */
	private int[] fitness;
 
	/**
	 * ��Ⱥ�и���������ۼƸ���
	 */
	private float[] Pi; 
	
	/**
	 * �������
	 */
	private float Pc;
	
	/**
	 * �������
	 */
	private float Pm;
	
	/**
	 * ��ǰ����
	 */
	private int t;
 
	private Random random;
 
	public GeneticAlgorithmTSP() {}
 
	/**
	 * �������в���
	 * @param scale ��Ⱥ��ģ
	 * @param cityNum ����������Ⱦɫ�峤��
	 * @param max_gen ��ֳ����
	 * @param Pc �������
	 * @param Pm �������
	 */
	public GeneticAlgorithmTSP(int scale, int cityNum, int max_gen, float Pc, float Pm) {
		this.scale = scale;
		this.cityNum = cityNum;
		MAX_GEN = max_gen;
		this.Pc = Pc;
		this.Pm = Pm;
	}
 
	/**
	 * ��ʼ��GA�㷨��
	 * @param filename �����ļ�·�������ļ��洢���г��нڵ���������
	 * @throws IOException
	 */
	public void init(String filename) throws IOException {
		// ��ȡ����
		int[] x = new int[cityNum];
		int[] y = new int[cityNum];
		distances = new int[cityNum][cityNum];
		BufferedReader data = new BufferedReader(new InputStreamReader(
                new FileInputStream(filename)));
		String strbuff = null;
		for (int i = 0; i < cityNum; i++) {
			// ��ȡһ�����ݣ����ݸ�ʽΪ 1 6734 1453
			strbuff = data.readLine();
			// �ַ��ָ�
			if(strbuff != null) {
			    String[] strcol = strbuff.split(" ");
			    // x����
			    x[i] = Integer.valueOf(strcol[1]);
			    // y����
			    y[i] = Integer.valueOf(strcol[2]);
			}
		}
		data.close();
		
		// ����������
		// ��Ծ������⣬������㷽��Ҳ��һ�����˴��õ���att48��Ϊ����������48�����У�������㷽��Ϊαŷ�Ͼ��룬����ֵΪ10628
		for (int i = 0; i < cityNum - 1; i++) {
			distances[i][i] = 0; // �Խ���Ϊ0���Լ����Լ��ľ���
			for (int j = i + 1; j < cityNum; j++) {
				double rij = Math
						.sqrt(((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j])
								* (y[i] - y[j])) / 10.0);
				// ��С������ȡ����һλ
				int tij = (int) Math.round(rij);// ��������ȡ��
				if (tij < rij) {
					distances[i][j] = tij + 1;
					distances[j][i] = distances[i][j];
				} else {
					distances[i][j] = tij;
					distances[j][i] = distances[i][j];
				}
			}
		}
		distances[cityNum - 1][cityNum - 1] = 0;
 
		bestLength = Integer.MAX_VALUE;
		bestTour = new int[cityNum + 1];
		bestT = 0;
		t = 0;
 
		newPopulation = new int[scale][cityNum];
		oldPopulation = new int[scale][cityNum];
		fitness = new int[scale];
		Pi = new float[scale];
 
		random = new Random(System.currentTimeMillis());
		
		/*
		 * for(int i=0;i<cityNum;i++) { for(int j=0;j<cityNum;j++) {
		 * System.out.print(distance[i][j]+","); } System.out.println(); }
		 */
	}
 
	/**
	 * ��ʼ����Ⱥ
	 */
	public void initPopulation() {
		int i, j, k;
		// Random random = new Random(System.currentTimeMillis());
		// ��ϻ���ȷ��ÿ�����壨Ⱦɫ�壩�ϵĻ���Ƭ�β��ظ�����֤Ѳ��·����ʵ�����壩
		for (k = 0; k < scale; k++)
		{
			oldPopulation[k][0] = random.nextInt(65535) % cityNum;
			for (i = 1; i < cityNum;)
			{
				oldPopulation[k][i] = random.nextInt(65535) % cityNum;
				for (j = 0; j < i; j++) {
					if (oldPopulation[k][i] == oldPopulation[k][j]) {
						break;
					}
				}
				if (j == i) {
					i++;
				}
			}
		}
 
		/*
		 * for(i=0;i<scale;i++) { for(j=0;j<cityNum;j++) {
		 * System.out.print(oldPopulation[i][j]+","); } System.out.println(); }
		 */
	}
 
	/**
	 * �����Ⱦɫ��ȫ����Ѳ��·��ȫ����
	 * @param chromosome ���壨Ⱦɫ�壩
	 * @return Ѳ��·��αŷ�Ͼ���ȫ��
	 */
	public int evaluate(int[] chromosome) {
		int len = 0;
		// Ⱦɫ�壬��ʼ����,����1,����2...����n
		for (int i = 1; i < cityNum; i++) {
			len += distances[chromosome[i - 1]][chromosome[i]];
		}
		// ����n����ʼ����
		len += distances[chromosome[cityNum - 1]][chromosome[0]];
		return len;
	}
 
	/**
	 * ������Ⱥ�и���������ۻ����ʣ�ǰ�����Ѿ�����������������Ӧ��fitness[max]����Ϊ����ѡ�����һ���֣�Pi[max]
	 */
	public void countRate() {
		int k;
		// ��Ӧ���ܺ�
		double sumFitness = 0;
		double[] tempf = new double[scale];
 
		for (k = 0; k < scale; k++) {
			tempf[k] = 10.0 / fitness[k];//
			sumFitness += tempf[k];
		}
 
		Pi[0] = (float) (tempf[0] / sumFitness);
		for (k = 1; k < scale; k++) {
			Pi[k] = (float) (tempf[k] / sumFitness + Pi[k - 1]);
		}
 
		/*
		 * for(k=0;k<scale;k++) { System.out.println(fitness[k]+" "+Pi[k]); }
		 */
	}
 
	
	/**
	 *��ѡĳ����Ⱥ����Ӧ����ߵĸ��壬ֱ�Ӹ��Ƶ��Ӵ�
	 *ǰ�����Ѿ�����������������Ӧ��Fitness[max]
	 */
	public void selectBestGh() {
		int k, i, maxid;
		int maxevaluation;
 
		maxid = 0;
		maxevaluation = fitness[0];
		for (k = 1; k < scale; k++) {
			if (maxevaluation > fitness[k]) {
				maxevaluation = fitness[k];
				maxid = k;
			}
		}
 
		if (bestLength > maxevaluation) {
			bestLength = maxevaluation;
			// ��õ�Ⱦɫ����ֵĴ���;
			bestT = t;
			for (i = 0; i < cityNum; i++) {
				bestTour[i] = oldPopulation[maxid][i];
			}
		}
 
		// System.out.println("���� " + t + " " + maxevaluation);
		// ��������Ⱥ����Ӧ����ߵ�Ⱦɫ��k���Ƶ�����Ⱥ�У����ڵ�һλ0
		copyGh(0, maxid);
	}
 
	
	/**
	 * ����Ⱦɫ��
	 * @param k ��Ⱦɫ������Ⱥ�е�λ��
	 * @param kk �ɵ�Ⱦɫ������Ⱥ�е�λ��
	 */
	public void copyGh(int k, int kk) {
		int i;
		for (i = 0; i < cityNum; i++) {
			newPopulation[k][i] = oldPopulation[kk][i];
		}
	}
 
	/**
	 * ѡ�����ӣ�����ѡ��������̶ĵķ�ʽ��
	 * ѡ��scale - 1�������Ŵ�
	 */
	public void select() {
		int k, i, selectId;
		float ran1;
		// Random random = new Random(System.currentTimeMillis());
		for (k = 1; k < scale; k++) {
		    // ����0��1֮����ұ�����λС���������
			ran1 = (float) (random.nextInt(65535) % 1000 / 1000.0);
			// System.out.println("����"+ran1);
			// ����ѡ�����
			for (i = 0; i < scale; i++) {
				if (ran1 <= Pi[i]) {
					break;
				}
			}
			selectId = i;
			// System.out.print("ѡ��" + selectId);
			copyGh(k, selectId);
		}
	}
 
	/**
	 * �����������ʽ��н���ͱ������
	 */
	public void evolution() {
		int k;
		// ��ѡĳ����Ⱥ����Ӧ����ߵĸ�����ӵ��Ŵ�������
		selectBestGh();
 
		// ����ѡ�������ѡscale-1����һ������
		select();
 
		// Random random = new Random(System.currentTimeMillis());
		float r;
 
		// ����ͱ������
		for (k = 0; k < scale; k = k + 2) {
			r = random.nextFloat();// /��������
			// System.out.println("������..." + r);
			if (r < Pc) {
				// System.out.println(k + "��" + k + 1 + "���н���...");
				//OXCross(k, k + 1);// ���н���
				OXCross1(k, k + 1);
			} else {
				r = random.nextFloat();// /��������
				// System.out.println("������1..." + r);
				// ����
				if (r < Pm) {
					// System.out.println(k + "����...");
					OnCVariation(k);
				}
				r = random.nextFloat();// /��������
				// System.out.println("������2..." + r);
				// ����
				if (r < Pm) {
					// System.out.println(k + 1 + "����...");
					OnCVariation(k + 1);
				}
			}
 
		}
	}
 
	/**
     * �����������ʽ��н���ͱ������
     * �������Ⱦɫ�岻���н������
     */
	public void evolution1() {
		int k;
		// ��ѡĳ����Ⱥ����Ӧ����ߵĸ�����ӵ��Ŵ�������
		selectBestGh();
 
		// ����ѡ�������ѡscale-1����һ������
		select();
 
		// Random random = new Random(System.currentTimeMillis());
		float r;
 
		for (k = 1; k + 1 < scale / 2; k = k + 2) {
			r = random.nextFloat();// /��������
			if (r < Pc) {
				OXCross1(k, k + 1);// ���н���
				//OXCross(k,k+1);//���н���
			} else {
				r = random.nextFloat();// /��������
				// ����
				if (r < Pm) {
					OnCVariation(k);
				}
				r = random.nextFloat();// /��������
				// ����
				if (r < Pm) {
					OnCVariation(k + 1);
				}
			}
		}
		if (k == scale / 2 - 1)// ʣ���һ��Ⱦɫ��û�н���L-1
		{
			r = random.nextFloat();// /��������
			if (r < Pm) {
				OnCVariation(k);
			}
		}
 
	}
 
	/**
	 * ��������
	 * @param k1 ����1
	 * @param k2 ����2
	 */
	public void OXCross(int k1, int k2) {
		int i, j, k, flag;
		int ran1, ran2, temp;
		int[] Gh1 = new int[cityNum];
		int[] Gh2 = new int[cityNum];
		// Random random = new Random(System.currentTimeMillis());
 
		ran1 = random.nextInt(65535) % cityNum;
		ran2 = random.nextInt(65535) % cityNum;
		// System.out.println();
		// System.out.println("-----------------------");
		// System.out.println("----"+ran1+"----"+ran2);
 
		while (ran1 == ran2) {
			ran2 = random.nextInt(65535) % cityNum;
		}
 
		if (ran1 > ran2)// ȷ��ran1<ran2
		{
			temp = ran1;
			ran1 = ran2;
			ran2 = temp;
		}
		// System.out.println();
		// System.out.println("-----------------------");
		// System.out.println("----"+ran1+"----"+ran2);
		// System.out.println("-----------------------");
		// System.out.println();
		flag = ran2 - ran1 + 1;// ɾ���ظ�����ǰȾɫ�峤��
		for (i = 0, j = ran1; i < flag; i++, j++) {
			Gh1[i] = newPopulation[k2][j];
			Gh2[i] = newPopulation[k1][j];
		}
		// �ѽ���ֵi=ran2-ran1������
 
		for (k = 0, j = flag; j < cityNum;)// Ⱦɫ�峤��
		{
			Gh1[j] = newPopulation[k1][k++];
			for (i = 0; i < flag; i++) {
				if (Gh1[i] == Gh1[j]) {
					break;
				}
			}
			if (i == flag) {
				j++;
			}
		}
 
		for (k = 0, j = flag; j < cityNum;)// Ⱦɫ�峤��
		{
			Gh2[j] = newPopulation[k2][k++];
			for (i = 0; i < flag; i++) {
				if (Gh2[i] == Gh2[j]) {
					break;
				}
			}
			if (i == flag) {
				j++;
			}
		}
 
		for (i = 0; i < cityNum; i++) {
			newPopulation[k1][i] = Gh1[i];// ������ϷŻ���Ⱥ
			newPopulation[k2][i] = Gh2[i];// ������ϷŻ���Ⱥ
		}
 
		// System.out.println("���н���--------------------------");
		// System.out.println(k1+"�����...");
		// for (i = 0; i < cityNum; i++) {
		// System.out.print(newPopulation[k1][i] + "-");
		// }
		// System.out.println();
		// System.out.println(k2+"�����...");
		// for (i = 0; i < cityNum; i++) {
		// System.out.print(newPopulation[k2][i] + "-");
		// }
		// System.out.println();
		// System.out.println("�������--------------------------");
	}
 
	/**
	 * ��������,��ͬȾɫ�彻�������ͬ�Ӵ�Ⱦɫ��
	 * @param k1 ����1
	 * @param k2 ����2
	 */
	public void OXCross1(int k1, int k2) {
		int i, j, k, flag;
		int ran1, ran2, temp;
		int[] Gh1 = new int[cityNum];
		int[] Gh2 = new int[cityNum];
		// Random random = new Random(System.currentTimeMillis());
 
		ran1 = random.nextInt(65535) % cityNum;
		ran2 = random.nextInt(65535) % cityNum;
		while (ran1 == ran2) {
			ran2 = random.nextInt(65535) % cityNum;
		}
 
		if (ran1 > ran2)// ȷ��ran1<ran2
		{
			temp = ran1;
			ran1 = ran2;
			ran2 = temp;
		}
 
		// ��Ⱦɫ��1�еĵ��������Ƶ�Ⱦɫ��2���ײ�
		for (i = 0, j = ran2; j < cityNum; i++, j++) {
			Gh2[i] = newPopulation[k1][j];
		}
 
		flag = i;// Ⱦɫ��2ԭ����ʼλ��
 
		for (k = 0, j = flag; j < cityNum;)// Ⱦɫ�峤��
		{
			Gh2[j] = newPopulation[k2][k++];
			for (i = 0; i < flag; i++) {
				if (Gh2[i] == Gh2[j]) {
					break;
				}
			}
			if (i == flag) {
				j++;
			}
		}
 
		flag = ran1;
		for (k = 0, j = 0; k < cityNum;)// Ⱦɫ�峤��
		{
			Gh1[j] = newPopulation[k1][k++];
			for (i = 0; i < flag; i++) {
				if (newPopulation[k2][i] == Gh1[j]) {
					break;
				}
			}
			if (i == flag) {
				j++;
			}
		}
 
		flag = cityNum - ran1;
 
		for (i = 0, j = flag; j < cityNum; j++, i++) {
			Gh1[j] = newPopulation[k2][i];
		}
 
		for (i = 0; i < cityNum; i++) {
			newPopulation[k1][i] = Gh1[i];// ������ϷŻ���Ⱥ
			newPopulation[k2][i] = Gh2[i];// ������ϷŻ���Ⱥ
		}
	}
 
	/**
	 * ��ζԻ���������
	 * @param k Ҫ����ĸ���
	 */
	public void OnCVariation(int k) {
		int ran1, ran2, temp;
		int count;// �Ի�����
 
		// Random random = new Random(System.currentTimeMillis());
		count = random.nextInt(65535) % cityNum;
 
		for (int i = 0; i < count; i++) {
 
			ran1 = random.nextInt(65535) % cityNum;
			ran2 = random.nextInt(65535) % cityNum;
			while (ran1 == ran2) {
				ran2 = random.nextInt(65535) % cityNum;
			}
			temp = newPopulation[k][ran1];
			newPopulation[k][ran1] = newPopulation[k][ran2];
			newPopulation[k][ran2] = temp;
		}
 
		/*
		 * for(i=0;i<L;i++) { printf("%d ",newGroup[k][i]); } printf("\n");
		 */
	}
 
	/**
	 * �����Ŵ��㷨���TSP
	 */
	public void solve() {
		int i;
		int k;
 
		// ��ʼ����Ⱥ
		initPopulation();
		// �����ʼ����Ⱥ��Ӧ�ȣ�Fitness[max]
		for (k = 0; k < scale; k++) {
			fitness[k] = evaluate(oldPopulation[k]);
			// System.out.println(fitness[k]);
		}
		// �����ʼ����Ⱥ�и���������ۻ����ʣ�Pi[max]
		countRate();
		//System.out.println("��ʼ��Ⱥ...");
		/*for (k = 0; k < scale; k++) {
			for (i = 0; i < cityNum; i++) {
				System.out.print(oldPopulation[k][i] + ",");
			}
			System.out.println();
			System.out.println("----" + fitness[k] + " " + Pi[k]);
		}*/
		
		for (t = 0; t < MAX_GEN; t++) {
			evolution1();
			//evolution();
			// ������ȺnewGroup���Ƶ�����ȺoldGroup�У�׼����һ������
			for (k = 0; k < scale; k++) {
				for (i = 0; i < cityNum; i++) {
					oldPopulation[k][i] = newPopulation[k][i];
				}
			}
			// ������Ⱥ��Ӧ��
			for (k = 0; k < scale; k++) {
				fitness[k] = evaluate(oldPopulation[k]);
			}
			// ������Ⱥ�и���������ۻ�����
			countRate();
		}
 
		/*System.out.println("�����Ⱥ...");
		for (k = 0; k < scale; k++) {
			for (i = 0; i < cityNum; i++) {
				System.out.print(oldPopulation[k][i] + ",");
			}
			System.out.println();
			System.out.println("---" + fitness[k] + " " + Pi[k]);
		}*/
 
		System.out.println("��ѳ��ȳ��ִ�����");
		System.out.println(bestT);
		System.out.println("��ѳ���");
		System.out.println(bestLength);
		System.out.println("���·����");
		for (i = 0; i < cityNum; i++) {
			System.out.print(bestTour[i] + ",");
		}
 
	}
 
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		System.out.println("Start....");
		long start = System.currentTimeMillis();
		GeneticAlgorithmTSP ga = new GeneticAlgorithmTSP(30, 48, 1000, 0.8f, 0.9f);
		ga.init("./src/cn/edu/jxufe/agl/data.txt");
		ga.solve();
		long end = System.currentTimeMillis();
		System.out.println("\n��ʱ��" + (end-start) + "ms");
	}
 
}