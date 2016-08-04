\name{adjust_columns}
\alias{adjust_columns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Add Columns of a Dataframe Based on Another Dataframe
}
\description{
Compare column names of one dataframe to column names of the template dataframe. After that, add missing columns and fill them with NA values. 
}
\usage{
adjust_columns(df_to_adjust = NULL, df_template = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_to_adjust}{
  dataframe to be adjusted.

}
  \item{df_template}{
  dataframe to be used as a template.
}
}
\author{
Ivan Grishagin
}