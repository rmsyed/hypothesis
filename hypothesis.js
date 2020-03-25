class htest {


	/*......................... Helper functions......................... */
	/* Helper function: compute the mean of an array */
	static avg(x) {
		return this.sum(x) / x.length;
	}


	/* Helper function: perform rounding to N digits */
	static round(n, digits) {
		return Math.round(n * (10**digits)) / (10**digits);
	}	

	/* Helper function: compute the sum of an array */
	static sum(x) {
		var s=0;
		for(let i=0; i<x.length; i++) {
	    	s += x[i];
	  }
	  return s
	}


	/* Helper function: remove NaNs in the "table" of arrays x and y */
	static remove_nans(x,y) {
		/* identify NaN or non-numeric rows for exclusion */
		var bad_rows = new Set();
		var tempx=[], tempy=[], i;
		for(i=0; i<x.length; i++) {
			if (!(isNaN(x[i]) || isNaN(y[i]))) {
				tempx.push(x[i]);
				tempy.push(y[i]);
			}
		}
		return [tempx, tempy];
	}



	/* Initialize the desired rounding of digits */
	static init() {
		this.sig_level = 0.05;
		this.sig_digits = 4;
	}





	/*......................... Distribution functions......................... */
	/* Source: 	https://www.math.ucla.edu/~tom/distributions/ */


	/* Compute LogGamma using Lanczos' approximation */
	static LogGamma(Z) {
		var S = 1 + 76.18009173/Z - 86.50532033/(Z+1) + 24.01409822/(Z+2) - 1.231739516/(Z+3) + .00120858003/(Z+4) - .00000536382/(Z+5);
		var LG = (Z-.5) * Math.log(Z+4.5) - (Z+4.5) + Math.log(S*2.50662827465);
		return LG
	}

	/* Compute the Incomplete Beta function */
	static  Betinc(X,A,B) {
		var A0=0, B0=1, A1=1, B1=1, M9=0, A2=0, C9;
		while (Math.abs((A1-A2)/A1)>.0000001) {
			A2=A1;
			C9=-(A+M9)*(A+B+M9)*X/(A+2*M9)/(A+2*M9+1);
			A0=A1+C9*A0;
			B0=B1+C9*B0;
			M9=M9+1;
			C9=M9*(B-M9)*X/(A+2*M9-1)/(A+2*M9);
			A1=A0+C9*A1;
			B1=B0+C9*B1;
			A0=A0/B1;
			B0=B0/B1;
			A1=A1/B1;
			B1=1;
		}
		return A1/A
	}


	/* Beta function CDF */
	static Betacdf(Z,A,B) {
	    var S, BT, Bcdf;
		S = A + B;
		BT = Math.exp(this.LogGamma(S) - this.LogGamma(B) - this.LogGamma(A) + A * Math.log(Z) + B * Math.log(1-Z));
		if (Z < (A+1) / (S+2)) {
			Bcdf = BT * this.Betinc(Z,A,B)
		} else {
			Bcdf = 1 - BT * this.Betinc(1-Z, B, A)
		}
		
		return Bcdf
	}


	/* Compute the Gamma CDF (for X >= A + 1) */
	static Gcf(X,A) {
        var A0 = 0;
        var B0 = 1;
        var A1 = 1;
        var B1 = X;
        var AOLD = 0;
        var N = 0;
        while (Math.abs((A1 - AOLD)/A1)>.00001) {
            AOLD = A1;
            N = N+1;
            A0 = A1+(N-A)*A0;
            B0 = B1+(N-A)*B0;
            A1 = X*A0+N*A1;
            B1 = X*B0+N*B1;
            A0 = A0/B1;
            B0 = B0/B1;
            A1 = A1/B1;
            B1 = 1;
        }
        var Prob=Math.exp(A*Math.log(X)-X-this.LogGamma(A))*A1;
        
        return 1-Prob
    }


    /* Compute the Gamma CDF (for X < A + 1) */
    static Gser(X,A) {
        var T9 = 1/A;
        var G = T9;
        var I = 1;
        while (T9 > G*.00001) {
            T9 = T9 * X / (A + I);
            G = G + T9;
            I = I + 1;
        }
        G = G * Math.exp(A * Math.log(X) - X - this.LogGamma(A));
        return G
    }

    /* Compute the Gamma CDF */
    static Gammacdf(x,a) {
        var GI;
        if (x<=0) {
            GI = 0
        } else if (x < a + 1) {
            GI = this.Gser(x,a)
        } else {
            GI = this.Gcf(x,a)
        }
        return GI
    }


	/* Compute the Chi-squared CDF */
    static cdf_chisq(x, df) {
        let Z = x
        let DF = df;
        let Chisqcdf = this.Gammacdf(Z/2, DF/2)

        Chisqcdf = Math.round(Chisqcdf * 100000) / 100000;
        return 1 - Chisqcdf;
    }

	/* Compute the F distribution CDF */
	static cdf_f(f_stat, f1, f2) {
	    var X = f_stat;
		var Z = X / (X + f2/f1);
		var Fcdf = this.Betacdf(Z, f1/2, f2/2);
		Fcdf = Math.round(Fcdf * 100000) / 100000;
	    return Fcdf;
	}	

	/* Compute the t-distribution CDF */
	static cdf_t(t_stat, df) {
		var X = t_stat;
		var tcdf, betacdf;

		var A = df/2;
		var S = A+.5;
		var Z = df / (df + X*X);
		var BT = Math.exp(this.LogGamma(S)-this.LogGamma(.5)-this.LogGamma(A) + A * Math.log(Z) + .5 * Math.log(1-Z));
		if (Z < (A+1) / (S+2)) {
			betacdf = BT * this.Betinc(Z,A,.5)
		} else {
			betacdf = 1 - BT * this.Betinc(1-Z,.5,A)
		}
		if (X<0) {
			tcdf = betacdf/2
		} else {
			tcdf = 1 - betacdf/2
		}
		tcdf = Math.round(tcdf*100000)/100000;
		
	    return tcdf;
	}




    /*......................... Hypothesis Testing functions......................... */



    /**
    * @typedef {Object} ChisqResult
    * @property {number} statistic - Chi-squared statistic
    * @property {number} df - degrees of freedom
    * @property {number} pval - p-value
    */

    /**
    * @typedef {Object} AnovaResult
    * @property {number} statistic - f-statistic
    * @property {number} df_1 - degrees of freedom (within)
    * @property {number} df_2 - degrees of freedom (between)
    * @property {number} pval - p-value
    */

    /**
	* @typedef {Object} CorrelationResult
	* @property {number} statistic - the test statistic
	* @property {number} df - degrees of freedom
	* @property {number} t_stat - t-statistic
	* @property {number} pval - p-value
	*/







    /**
	* Run the Chi-squared test
	* @param {array} x - Array of length N (categorical values)
    * @param {array} y - Array of length N (categorical values)
    * @returns {ChisqResult}
    */

    static chisq(x,y) {
    	/* Assumes both x and y are categoricals */
    	/* Maintain a table of the counts of each. So for each unique y value
    	what is the count of each unique x value.*/

        var i, k;
        let counts = {};
        var allkeys = new Set();
        for(i=0; i<y.length; i++) {
        	allkeys.add(x[i]);
        	if(!(y[i] in counts)) {
        		counts[y[i]] = {};
        	}

    		if(!(x[i] in counts[y[i]])) {
    			counts[y[i]][x[i]] = 1
    		}
    		else {
    			counts[y[i]][x[i]] += 1;
    		}
        }


        /* inject zeros for missing keys */
        allkeys = Array.from(allkeys);
        for(i=0; i<allkeys.length; i++) {
        	for(const y_val in counts) {
        		if(!(allkeys[i] in counts[y_val])) {
        			counts[y_val][allkeys[i]] = 0;
        		}
        	}
        }

        /* now flatten this dictionary into a 2D array; */
        var frame = [];
        for(const y_val in counts) {
        	var temp = [];
        	for(i=0; i<allkeys.length; i++) {
        		temp.push(counts[y_val][allkeys[i]]);
        	}
        	frame.push(temp);
        }



        /* compute the columnwise and row-wise aggregates */
        let rowaggs = {};
        let colaggs = {};
        let total = 0;

        for(i=0; i<frame.length; i++) {
            rowaggs[i] = 0;
            for(k=0; k<frame[i].length; k++) {
                rowaggs[i] += frame[i][k];
                total += frame[i][k];
            }
        }

        for(i=0; i<frame[0].length; i++) {
            colaggs[i] = 0;
            for(k=0; k<frame.length; k++) {
                colaggs[i] += frame[k][i];
            }
        }

        /* for each cell, compute what the total fraction would be for this column */
        let chi_error = 0;
        for(i=0; i<frame.length; i++) {
            let rowfrac = rowaggs[i] / total;
            for(k=0; k<frame[i].length; k++) {
                let expects = rowfrac * colaggs[k];
                chi_error += Math.pow(expects - frame[i][k], 2) / expects;
            }
        }        

        /* Return the results */
        var results = this.cdf_chisq(chi_error, (frame.length - 1) * (frame[0].length - 1));
        return {
	        "statistic": this.round(chi_error, this.sig_digits),
	        "df": this.round((frame.length - 1) * (frame[0].length - 1), this.sig_digits),
	        "pval": this.round(results, this.sig_digits)
        }
    }



    /**
	* Run the one-way ANOVA test
	* @param {array} x - Array of length N (continuous values)
    * @param {array} y - Array of length N (categorical values)
    * @returns {AnovaResult}
    */

	static anova(x, y) {

		/* Find the unique categories in y */
		var cats = new Set();
		var i,k,m, df1, df2;
		for(i=0; i<y.length; i++) {
			cats.add(y[i]);
		}
		cats = Array.from(cats);
		var cat_slivers = {};
		for(i=0; i<cats.length; i++) {
			cat_slivers[cats[i]] = [];
		}

		/* separate x into category arrays */
		for(i=0; i<x.length; i++) {
			cat_slivers[y[i]].push(x[i]);
		}

		/* calculate the within-category averages */
		var cat_avgs = {};
		for(i=0; i<cats.length; i++) {
			let cat = cats[i];
			var avg = 0;
			cat_avgs[cat] = this.avg(cat_slivers[cat]);
		}

		/* Calculate the weighted averages of the between-group averages */
		var cat_avgs_global = 0;
		for(i=0; i<cats.length; i++) {
			let cat = cats[i];
			cat_avgs_global += cat_avgs[cat] * cat_slivers[cat].length;
		}
		cat_avgs_global = cat_avgs_global / y.length;

		/* calculate the between-groups sum of squared differences */
		var SB = 0;
		for(i=0; i<cats.length; i++) {
			let cat = cats[i];
			SB += cat_slivers[cat].length * (cat_avgs[cat] - cat_avgs_global)**2;
		}
		df1 = cats.length - 1;
		var MSB = SB / df1 // mean-squared between-groups value




		/* calculate the within-groups sum of squared differences */
		var SW = 0;
		for(i=0; i<cats.length; i++) {
			let cat = cats[i];
			var curarr = cat_slivers[cat];
			for(k=0; k<curarr.length; k++) {
				SW += (curarr[k] - cat_avgs[cat])**2
			}
		}
		df2 = 0;
		for(i=0; i<cats.length; i++) {
			let cat = cats[i];
			df2 += cat_slivers[cat].length - 1
		}
		var MSW = SW / df2 // mean-squared within-groups values


		/* calculate the F-statistic */
		var f_stat = MSB / MSW;

		let cdf_val = this.cdf_f(f_stat, df1, df2);

		/* get the corresponding p-value now. */

		cdf_val = (1 - cdf_val);
		return {
			"statistic": this.round(f_stat, this.sig_digits),
			"df_1": this.round(df1, this.sig_digits),
			"df_2": this.round(df2, this.sig_digits),
			"pval": this.round(cdf_val, this.sig_digits)
		};
	}




	/**
	* Run the Spearman rank correlation test
	* @param {array} x - Array of length N (continuous values)
    * @param {array} y - Array of length N (continuous values)
    * @returns {CorrelationResult}
    */
	static spearman(x, y) {
		var tuples = [];
		var i,k,m;

		/* get rid of NaN values */
		var res = this.remove_nans(x,y);
		x = res[0];
		y = res[1];

		/* sort the values to get the ranks */
		var tx = x.slice(0);
		tx.sort(function(a,b) {return parseFloat(a) - parseFloat(b)}) // got to specify this to make sure we are sorting numerically
		var ty = y.slice(0);
		ty.sort(function(a,b) {return parseFloat(a) - parseFloat(b)});

		/* compute the ranks */
		// for each sorted number, we will assign it a rank. If the rank already exists, we will add the new rank value
		// and average these at the end.
		var rank_map = {};
		var arrs = [tx, ty];
		var has_ties = false;

		for(k=0; k<arrs.length; k++) {
			rank_map = {};
			var curarr = arrs[k];
			for(i=0; i<curarr.length; i++) {
				if(curarr[i] in rank_map) {
					rank_map[curarr[i]]["total"] += (i+1);
					rank_map[curarr[i]]["count"] += 1;
					has_ties = true;
				}
				else {
					rank_map[curarr[i]] = {"total":(i+1), "count":1};
				}
			}

			/* average the tied ranks */
			for (const prop in rank_map) {
				rank_map[prop] = rank_map[prop]["total"] / rank_map[prop]["count"];
			}

			/* map the raw input values to their corresponding ranks */
			for(i=0; i<curarr.length; i++) {
				if(k==0) {
					tx[i] = rank_map[x[i]];
				}
				else {
					ty[i] = rank_map[y[i]];
				}
			}
		}

		/* calculate the S-statistic */
		var s_stat = 0;
		for(i=0; i<tx.length; i++) {
			s_stat += (tx[i] - ty[i])**2;
		}
		s_stat = 1 - ((s_stat * 6) / (tx.length * (tx.length**2 - 1)))



		/* if there are ties, go with pearson. */
		var results = [];
		if(has_ties) {
			results = this.pearson(tx, ty);
		}
		else {
			/* compute the p-value */
			let t_stat = Math.sqrt((x.length - 2) / (1 - s_stat**2)) * s_stat;
			let cdf_val = this.cdf_t(t_stat, x.length - 2);
			if(t_stat < 0) {
				cdf_val = cdf_val * 2;
			}
			else {
				cdf_val = (1 - cdf_val)*2;
			}
			//results = [s_stat, t_stat, cdf_val];
			results = {
				"statistic": this.round(s_stat, this.sig_digits),
				"df": this.round(x.length - 2, this.sig_digits),
				"pval": this.round(cdf_val, this.sig_digits),
				"t_stat": this.round(t_stat, this.sig_digits)
			}
		}
		return results
	}



	/**
	* Run the Pearson's correlation test
	* @param {array} x - Array of length N (continuous values)
    * @param {array} y - Array of length N (continuous values)
    * @returns {CorrelationResult}
    */
	static pearson(x, y) {
		//this.init();
		if(x.length != y.length) {
			return [null, null];
		}
		var res = this.remove_nans(x,y);
		x = res[0];
		y = res[1];

		let i=0;
		let numer = 0, denom = 0;
		let x_avg = x.reduce((prevval, curval) => prevval + curval, 0) / x.length;
		let y_avg = y.reduce((prevval, curval) => prevval + curval, 0) / y.length;

		for(i=0; i<x.length; i++) {
			x[i] = x[i] - x_avg;
			y[i] = y[i] - y_avg;
		}
		for(i=0; i<x.length; i++) {
			numer += x[i] * y[i];
		}
		for(i=0; i<x.length; i++) {
			denom += x[i]**2
		}
		denom = Math.sqrt(denom)
		let d1 = denom;
		denom = 0;
		for(i=0; i<x.length; i++) {
			denom += y[i]**2
		}
		denom = d1 * Math.sqrt(denom);

		let rho = numer / denom;
		let t_stat = Math.sqrt((x.length - 2) / (1 - rho**2)) * rho;

		// need to get the corresponding p-value now.
		let cdf_val = this.cdf_t(t_stat, x.length - 2);
		if(t_stat < 0) {
			cdf_val = cdf_val * 2;
		}
		else {
			cdf_val = (1 - cdf_val)*2;
		}
		return {
			"statistic": this.round(rho, this.sig_digits),
			"df": this.round(x.length - 2, this.sig_digits),
			"pval": this.round(cdf_val, this.sig_digits),
			"t_stat": this.round(t_stat, this.sig_digits)
			}
	}
}

stats.init()
