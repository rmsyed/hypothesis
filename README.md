# hypothesis
Easy-to-use javascript library for statistical hypothesis testing. Your code only needs to pass the arrays on which to perform the tests and **hypothesis** does the rest!

```html
<script type="text/javascript src="hypothesis.js"></script>
```

# Tests supported
#### Currently, hypothesis supports the following tests:
- Correlation tests (Pearson and Spearman)
- Chi-squared test
- One-way ANOVA test

# Usage
### Pearson correlation test
```javascript
var x = [4,4,1,1,2,1,4,2,2,4,4,3,3,3,4,4,4,1,2,1,1,2,2,4,2,1,2,2,4,6,8,2];
var y = [160,160,108,258,360,225,360,146.7,140.8,167.6,167.6,275.8,275.8,275.8,472,460,440,78.7,75.7,71.1,120.1,318,304,350,400,79,120.3,95.1,351,145,301,121];
stats.pearson(x, y);
```
##### Returns:
```javascript
{
  "df": 30,
  "pval": 0.0253,
  "statistic": 0.395,
  "t_stat": 2.3548
}
```


### Spearman correlation test
```javascript
var x = [4,4,1,1,2,1,4,2,2,4,4,3,3,3,4,4,4,1,2,1,1,2,2,4,2,1,2,2,4,6,8,2];
var y = [160,160,108,258,360,225,360,146.7,140.8,167.6,167.6,275.8,275.8,275.8,472,460,440,78.7,75.7,71.1,120.1,318,304,350,400,79,120.3,95.1,351,145,301,121];
stats.spearman(x, y);
```
##### Returns:
```javascript
{
  "df": 30,
  "pval": 0.0014,
  "statistic": 0.5398,
  "t_stat": 3.5121
}
```


### One-way ANOVA test
```javascript
var x = [160,160,108,258,360,225,360,146.7,140.8,167.6,167.6,275.8,275.8,275.8,472,460,440,78.7,75.7,71.1,120.1,318,304,350,400,79,120.3,95.1,351,145,301,121];
var y = [4,4,1,1,2,1,4,2,2,4,4,3,3,3,4,4,4,1,2,1,1,2,2,4,2,1,2,2,4,6,8,2];
stats.anova(x, y); // y must be the categorical array
```
##### Returns:
```javascript
{
  "df_1": 5,
  "df_2": 26,
  "statistic": 2.3817,
  "pval": 0.0662
}
```


### Chi-squared test
```javascript
var x = [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1];
var y = [4,4,1,1,2,1,4,2,2,4,4,3,3,3,4,4,4,1,2,1,1,2,2,4,2,1,2,2,4,6,8,2];
stats.chisq(x, y); // both x and y must be categorical arrays.
```
##### Returns:
```javascript
{
  "df": 5,
  "statistic": 6.2371,
  "pval": 0.2838
}
```
